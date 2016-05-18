#include<iostream>
#include<vector>
#include<fstream>
#include<string.h>
#include <GL/freeglut.h>
#include <GL/glui.h>
#include "Physika_Core/Vectors/vector.h"
#include "Physika_Core/Utilities/physika_assert.h"
#include "Physika_Render/OpenGL_Primitives/opengl_primitives.h"
#include "Physika_GUI/Glut_Window/glut_window.h"
#include "Physika_GUI/Glui_Window/glui_window.h"
#include "Physika_IO/Volumetric_Mesh_IO/volumetric_mesh_io.h"
#include "Physika_IO/Surface_Mesh_IO/surface_mesh_io.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/cubic_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/quad_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/tet_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/tri_mesh.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh_interpolation.h"
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"
#include "Physika_Render/Volumetric_Mesh_Render/volumetric_mesh_render.h"
#include "Physika_Core/Transform/transform.h"
#include "Physika_Dynamics/FEM/FEM_Solid_Force_Model/fem_solid_force_model.h"
#include "Physika_Dynamics/FEM/FEM_Solid_Force_Model/tri_tet_mesh_fem_solid_force_model.h"
#include "Physika_Dynamics/Constitutive_Models/st_venant_kirchhoff.h"
#include "Physika_Numerics/Linear_System_Solvers/linear_system.h"
#include "Physika_Numerics/Linear_System_Solvers/conjugate_gradient_solver.h"
#include "Physika_Numerics/Linear_System_Solvers/plain_generalized_vector_T.h"


using namespace std;
using namespace Physika;

template <typename Scalar, int Dim>
class StVKStiffness : public LinearSystem<Scalar>{
public:
	StVKStiffness(TriTetMeshFEMSolidForceModel<Scalar, Dim> &forcemodel, vector<unsigned int> &fixedPoints, vector<Vector<double, 3>> &cur_pos){
		forcemodel_ = &forcemodel;
		fixedPoints_ = &fixedPoints;
		cur_pos_ = &cur_pos;
	}
	virtual ~StVKStiffness(){

	}
	virtual void multiply(const GeneralizedVector<Scalar> &x, GeneralizedVector<Scalar> &result) const{
		vector<Vector<Scalar, Dim>> *dis_differential = new vector<Vector<Scalar, Dim>>(x.size() / 3);
		vector<Vector<Scalar, Dim>> *force_differential = new vector<Vector<Scalar, Dim>>(x.size() / 3);
		this->GeneralizedVec2vec3D(x, *dis_differential);
		forcemodel_->computeGlobalInternalForceDifferentials(*cur_pos_, *dis_differential, *force_differential);
		this->vec3D2GeneralizedVec(*force_differential, result);
		delete dis_differential;
		delete force_differential;
	}
	virtual void preconditionerMultiply(const GeneralizedVector<Scalar> &x, GeneralizedVector<Scalar> &result) const{

	}
	virtual double innerProduct(const GeneralizedVector<Scalar> &x, const GeneralizedVector<Scalar> &y) const{
		PHYSIKA_ASSERT(x.size() == y.size());
		const PlainGeneralizedVector<Scalar> &plain_x = dynamic_cast<const PlainGeneralizedVector<Scalar>&>(x);
		const PlainGeneralizedVector<Scalar> &plain_y = dynamic_cast<const PlainGeneralizedVector<Scalar>&>(y);
		double sum = 0;
		for (int i = 0; i < x.size(); ++i)sum += plain_x[i] * plain_y[i];
		return sum;
	}
	virtual void filter(GeneralizedVector<Scalar> &x) const {
		PlainGeneralizedVector<Scalar> &plain_x = dynamic_cast<PlainGeneralizedVector<Scalar>&> (x);
		for (int i = 0; i < (fixedPoints_->size()); ++i){
			plain_x[((*fixedPoints_)[i]) * 3] = 0;
			plain_x[((*fixedPoints_)[i]) * 3 + 1] = 0;
			plain_x[((*fixedPoints_)[i]) * 3 + 1] = 0;
		}
	}
	vector<Vector<double, 3>> * cur_pos_;
	TriTetMeshFEMSolidForceModel<Scalar, Dim> * forcemodel_;
	vector<unsigned int> *fixedPoints_;
protected:
	void vec3D2GeneralizedVec(const vector<Vector<Scalar, Dim>> &a, GeneralizedVector<Scalar>&b)const {
		PlainGeneralizedVector<Scalar> &plain_b = dynamic_cast<PlainGeneralizedVector<Scalar>&>(b);
		for (int i = 0; i < a.size(); ++i){
			plain_b[i * 3] = a[i][0];
			plain_b[i * 3 + 1] = a[i][1];
			plain_b[i * 3 + 2] = a[i][2];
		}
	}
	void GeneralizedVec2vec3D(const GeneralizedVector<Scalar>&b, vector<Vector<Scalar, Dim>> &a)const{
		const PlainGeneralizedVector<Scalar> &plain_b = dynamic_cast<const PlainGeneralizedVector<Scalar>&>(b);
		for (int i = 0; i < a.size(); ++i)a[i] = Vector<Scalar, Dim>(plain_b[3 * i], plain_b[3 * i + 1], plain_b[3 * i + 2]);
	}
};

vector<unsigned int> fixPoints;
PlainGeneralizedVector<double> force;
PlainGeneralizedVector<double> df;
PlainGeneralizedVector<double> dx;
VolumetricMeshRender<double, 3> meshRender;
VolumetricMeshRender<double, 3> meshRender2;
VolumetricMesh<double, 3> *vMesh;
vector<Vector<double, 3>> cur_pos;
vector<Vector<double, 3>> displacement;
StVKStiffness<double, 3> * plinearSys = NULL;
ConjugateGradientSolver<double> * psolver = NULL;
TriTetMeshFEMSolidForceModel<double, 3> *pforceModel = NULL;
GlutWindow glut_window;

void displayFunction()
{
	cout << "display" << endl;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);		     // Clear Screen and Depth Buffer
	GlutWindow *cur_window = (GlutWindow*)glutGetWindowData();
	//cur_window->orbitCameraRight(0.1);
	(cur_window->camera()).look();

	meshRender.disableRenderSolid();
	meshRender2.disableRenderSolid();
	//meshRender.disableRenderWireframe();
	//meshRender.disableRenderVertices();
	//meshRender.enableRenderVertices();
	meshRender.enableRenderWireframe();
	meshRender2.enableRenderWireframe();
	//meshRender.enableRenderSolid();
	//meshRender2.enableRenderSolid();	
	meshRender.setVolumetricMesh(vMesh, &displacement);
	meshRender2.setVolumetricMesh(vMesh);

	//meshRender2.renderVertexWithColor<double>(fixPoints, Color<double>(1, 0, 0));
	meshRender.renderVertexWithColor<double>(fixPoints, Color<double>(1, 0, 0));
	glColor4f(0, 1, 0, 1.0);
	meshRender2.render();
	glColor4f(1, 1, 1, 1.0);
	meshRender.render();
	//meshRender.renderVertexWithColor(vertex_id,color);
	//meshRender.renderElementWithColor(element_id,color);
	//glTranslatef(3,0,0);
	//meshRender2.renderSolidWithAlpha(0.4);

	cur_window->displayFrameRate();
	glutSwapBuffers();
	glFlush();
}

void idleFunction()
{/*
	vector<Vector<double, 3>> cur_force;
	pforceModel->computeGlobalInternalForces(cur_pos, cur_force);
	for (unsigned int i = 0; i < cur_force.size(); ++i){
		df[i * 3] = force[i * 3] - cur_force[i][0];
		df[i * 3 + 1] = force[i * 3 + 1] - cur_force[i][1];
		df[i * 3 + 2] = force[i * 3 + 2] - cur_force[i][2];
	}
	for (unsigned int i = 0; i < dx.size(); ++i)dx[i] = 0;
	psolver->solve(*plinearSys, df, dx);
	cout << "PCG Iterations:" << psolver->iterationsUsed() << endl;
	//cout << dx << endl;
	for (int i = 0; i < cur_pos.size(); ++i){
		cur_pos[i][0] += dx[i * 3];
		cur_pos[i][1] += dx[i * 3 + 1];
		cur_pos[i][2] += dx[i * 3 + 2];
		displacement[i] = cur_pos[i] - vMesh->vertPos(i);
	}
	psolver->reset();
	cout << "iteration :" << endl;*/
}

void initFunction()
{
	glClearColor(0.0, 0.0, 0.0, 1.0);
	glClearDepth(1.0);


}

void keyboardFunction(unsigned char key, int x, int y)
{
	GlutWindow::bindDefaultKeys(key, x, y);
	switch (key)
	{
	case 't':
		cout << "test\n";
		break;
	case 's':
		cout << "save picture!" << endl;
		glut_window.saveScreen();
		break;
	default:
		break;
	}
}

int main(){
	vMesh = VolumetricMeshIO<double, 3>::load(string("FEMTest/bar.smesh"));    //mesh information
	SurfaceMesh<double> sMesh;
	SurfaceMeshIO<double>::load(string("FEMTest/bar-render.obj"),&sMesh);
	VolumetricMeshInterpolation<double, 3> interpolation(*vMesh);
	interpolation.getSurfaceMeshWeights(sMesh);


	fstream filein("FEMTest/fixedPoints.txt");    //fixed points
	int num;
	cout << "fixed Points:" << endl;
	while (filein >> num){
		fixPoints.push_back(num - 1);
		cout << num - 1 << ' ';
	}
	cout << endl;
	filein.close();

	filein.open("FEMTest/constantForce.txt");   //constant force
	force.resize(vMesh->vertNum() * 3);
	dx.resize(vMesh->vertNum() * 3);
	for (unsigned int i = 0; i < vMesh->vertNum() * 3; ++i){ filein >> force[i]; force[i] = -force[i]; }
	filein.close();

	StVK<double, 3> stvk(1e6, 0.3, IsotropicHyperelasticMaterialInternal::YOUNG_AND_POISSON);
	vector<ConstitutiveModel<double, 3>*> constitutiveModels;
	constitutiveModels.push_back(&stvk);
	TriTetMeshFEMSolidForceModel<double, 3> forceModel(*vMesh, constitutiveModels);

	//计算每个顶点的初始位置向量   初始化
	for (unsigned int i = 0; i < vMesh->vertNum(); ++i){
		cur_pos.push_back(vMesh->vertPos(i));
		displacement.push_back(Vector<double, 3>(0, 0, 0));
		dx[i * 3] = 0;
		dx[i * 3 + 1] = 0;
		dx[i * 3 + 2] = 0;
	}
	df = force;

	StVKStiffness<double, 3> linearSys(forceModel, fixPoints, cur_pos);
	ConjugateGradientSolver<double> solver;
	plinearSys = &linearSys;
	psolver = &solver;
	pforceModel = &forceModel;
	solver.enableStatusLog();

	//render project
	glut_window.setCameraPosition(Vector<double, 3>(0, -5, 5));
	glut_window.setCameraFocusPosition(Vector<double, 3>(0, 0, 0));
	glut_window.setCameraNearClip(0.1);
	glut_window.setCameraFarClip(1.0e4);
	glut_window.setDisplayFunction(displayFunction);
	glut_window.setInitFunction(initFunction);
	cout << "Test GlutWindow with custom display function:\n";
	//glut_window.setIdleFunction(idleFunction);
	//glut_window.createWindow();
	//glut_window.mainLoopEvent();
	//glutPostRedisplay();
	//getchar();
	int j = 0;
	do
	{
		//glut_window.mainLoopEvent();
		//glutPostRedisplay();
		vector<Vector<double, 3>> cur_force;
		forceModel.computeGlobalInternalForces(cur_pos, cur_force);
		for (unsigned int i = 0; i < cur_force.size(); ++i){
			df[i * 3] = force[i * 3] - cur_force[i][0];
			df[i * 3 + 1] = force[i * 3 + 1] - cur_force[i][1];
			df[i * 3 + 2] = force[i * 3 + 2] - cur_force[i][2];
		}
		linearSys.filter(df);
		double dfsq = linearSys.innerProduct(df, df);
		cout << "df*df:" << dfsq << endl;
		if (dfsq < 1e-6) break;
		for (unsigned int i = 0; i < dx.size(); ++i)dx[i] = 0;
		solver.solve(linearSys, df, dx);
		cout << "PCG Iterations:" << solver.iterationsUsed() << endl;
		//cout << dx << endl;
		for (int i = 0; i < cur_pos.size(); ++i){
			cur_pos[i][0] += dx[i * 3];
			cur_pos[i][1] += dx[i * 3 + 1];
			cur_pos[i][2] += dx[i * 3 + 2];
			displacement[i] = cur_pos[i] - vMesh->vertPos(i);
		}
		solver.reset();
		cout << "iteration :" << ++j << endl;
		
	} while (j < 100);

	
	interpolation.interpolate(cur_pos, sMesh);
	SurfaceMeshIO<double>::save(string("FEMTest/bar_fine_d.obj"), &sMesh);
	cout << "save file OK" << endl;

	delete vMesh;
	return 0;
}

void mean_error(string file1, string file2){
	SurfaceMesh<double> sMesh1,sMesh2;
	SurfaceMeshIO<double>::load(string("FEMTest/bar-render.obj"), &sMesh1);
	SurfaceMeshIO<double>::load(string("FEMTest/bar_d.obj"), &sMesh2);
}