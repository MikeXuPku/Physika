#include<iostream>
#include "Physika_Geometry/Boundary_Meshes/surface_mesh.h"
#include "Physika_IO/Volumetric_Mesh_IO/volumetric_mesh_io.h"
#include "Physika_IO/Surface_Mesh_IO/surface_mesh_io.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh_interpolation.h"
#include "Physika_Dynamics/FEM/FEM_Solid_Force_Model/tri_tet_mesh_fem_solid_force_model.h"
#include "Physika_Dynamics/Constitutive_Models/st_venant_kirchhoff.h"
#include "FEMelastostatic.h"

VolumetricMesh<double, 3>* vMesh;
SurfaceMesh<double> sMesh;
vector<unsigned int> fixPoints;
PlainGeneralizedVector<double> force;
PlainGeneralizedVector<double> dx;
VolumetricMeshInterpolation<double, 3>* p_interpolation;
string filename;
using namespace std;
using namespace Physika;

void initialize(){
	filename = "elastostatic/bar-coarse";
	vMesh = VolumetricMeshIO<double, 3>::load(filename + string(".smesh"));    //mesh information
	SurfaceMeshIO<double>::load(string(filename + ".obj"), &sMesh);
	p_interpolation = new VolumetricMeshInterpolation<double, 3>(*vMesh);

	//p_interpolation->getSurfaceMeshWeights(sMesh);
	//p_interpolation->save(filename + "interpolation.txt");
	p_interpolation->load(filename+"interpolation.txt");
	fstream filein(filename + "Fixed.txt");    //fixed points
	int num;
	cout << "fixed Points:" << endl;
	while (filein >> num){
		fixPoints.push_back(num - 1);
		cout << num - 1 << ' ';
	}
	cout << endl;
	filein.close();

	filein.open(filename + "Force.txt");   //constant force
	force.resize(vMesh->vertNum() * 3);
	for (unsigned int i = 0; i < vMesh->vertNum() * 3; ++i){ filein >> force[i]; force[i] = -force[i]; }
	filein.close();
}

int main(){
	cout << "test begin:" << endl;
	initialize();   //初始化操作包括: 完成 体网格 表面网格 固定点 固定外力的读取
	cout << "initialize down!" << endl;

	//NeoHookeanAnisoTerm<double, 3> neoHookeanAnisoTerm(1e6, 0.3, IsotropicHyperelasticMaterialInternal::YOUNG_AND_POISSON, anisoDirection, 0);
	StVK<double, 3> stvk(1e6, 0.3, IsotropicHyperelasticMaterialInternal::YOUNG_AND_POISSON);
	vector<ConstitutiveModel<double, 3>*> constitutiveModels;
	constitutiveModels.push_back(&stvk);
	TriTetMeshFEMSolidForceModel<double, 3> forceModel(*vMesh, constitutiveModels);

	FEMelastostatic<double, 3> static_driver(vMesh,&forceModel,&fixPoints,&force);
	vector<Vector<double, 3>> static_position;
	static_driver.getStaticPosition(static_position);

	cout << "simulation down!" << endl;
	p_interpolation->interpolate(static_position, sMesh);
	SurfaceMeshIO<double>::save(filename + string("_static.obj"), &sMesh);
	cout << "save file OK" << endl;

	return 0;

}