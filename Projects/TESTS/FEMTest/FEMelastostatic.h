#include<iostream>
#include<vector>
#include<string.h>
#include "Physika_Core/Vectors/vector.h"
#include "Physika_Core/Utilities/physika_assert.h"
#include "Physika_Numerics/Linear_System_Solvers/conjugate_gradient_solver.h"
#include "Physika_Numerics/Linear_System_Solvers/plain_generalized_vector_T.h"
#include "Physika_Dynamics/FEM/FEM_Solid_Force_Model/fem_solid_force_model.h"
#include "Physika_Geometry/Volumetric_Meshes/volumetric_mesh.h"
#include "ForceModelLinearSys.h"

using namespace std;
using namespace Physika;


template <typename Scalar, int Dim>
class FEMelastostatic{

public:
	FEMelastostatic(FEMSolidForceModel<Scalar, Dim>*, vector<unsigned int>*, PlainGeneralizedVector<Scalar>*);
	~FEMelastostatic();
	void setFixedPoints(vector<unsigned int>*);
	void setOutForce(PlainGeneralizedVector<Scalar>*);
	void getStaticPosition(vector<Vector<Scalar, Dim>> &);



private:
	vector<Vector<double, Dim>> cur_pos_;
	ForceModelLinearSys<Scalar, Dim>* linearSys_;
	vector<unsigned int>* fixPoints_;
	PlainGeneralizedVector<Scalar>* force_;
	VolumetricMesh<Scalar, Dim> *vMesh_;
	PlainGeneralizedVector<Scalar> dx_;
};

template <typename Scalar, int Dim>
FEMelastostatic<Scalar,Dim>::FEMelastostatic(FEMSolidForceModel<Scalar, Dim> *forceModel, vector<unsigned int>*fixPoints, PlainGeneralizedVector<Scalar>*force):
cur_pos_(),linearSys_(forceModel,fixPoints,&cur_pos_),fixPoints_(fixPoints),force_(force){
	vMesh_ = &(forceModel->simulation_mesh_);
	for (unsigned int i = 0; i < vMesh_->vertNum(); ++i){
		cur_pos_.push_back(vMesh->vertPos(i));  //对cur_pos进行初始化
		displacement.push_back(Vector<Scalar, Dim>(0));
		for (unsigned int j = 0; j < Dim; ++j){
			dx[i * Dim + j] = 0;
		}
	}
	df = force;     //初始状态下 df 就等于 合外力
}

template <typename Scalar, int Dim>
FEMelastostatic<Scalar, Dim>::~FEMelastostatic(){

}

template <typename Scalar, int Dim>
void FEMelastostatic<Scalar, Dim>::setFixedPoints(vector<unsigned int>* fixPoints){
	fixPoints_ = fixPoints
}

template <typename Scalar, int Dim>
void FEMelastostatic<Scalar, Dim>::setOutForce(PlainGeneralizedVector<Scalar>* force){
	force_ = force;
}

template <typename Scalar, int Dim>
void FEMelastostatic<Scalar, Dim>::getStaticPosition(vector<Vector<Scalar, Dim>> &staticPosition){

	ConjugateGradientSolver<Scalar> solver;
	solver.enableStatusLog();

	int j = 0;
	do
	{
		vector<Vector<Scalar, Dim>> cur_force;
		forceModel.computeGlobalInternalForces(cur_pos_, cur_force);
		for (unsigned int i = 0; i < cur_force.size(); ++i){
			for (unsigned int kk = 0; kk < Dim; ++kk){
				df[i*Dim + kk] = force[i*Dim + kk] - cur_force[i][kk];
			}
		}
		linearSys_->filter(df);
		double dfsq = linearSys_->innerProduct(df, df);
		cout << "df*df:" << dfsq << endl;
		if (dfsq < 1e-6) break;
		for (unsigned int i = 0; i < dx_.size(); ++i)dx_[i] = 0;
		solver.solve(linearSys_, df, dx_);
		cout << "PCG Iterations:" << solver.iterationsUsed() << endl;
		//cout << dx_ << endl;
		for (int i = 0; i < cur_pos_.size(); ++i){
			for (unsigned int kk = 0; kk < Dim; ++kk){
				cur_pos_[i][kk] += dx[i*Dim + kk];
			}
		}
		solver.reset();
		cout << "iteration :" << ++j << endl;
	} while (j < 100);

	staticPostion = cur_pos_;
}