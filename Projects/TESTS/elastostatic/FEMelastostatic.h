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
	FEMelastostatic(VolumetricMesh<Scalar,Dim>* ,FEMSolidForceModel<Scalar, Dim>*, vector<unsigned int>*, PlainGeneralizedVector<Scalar>*);
	~FEMelastostatic();
	void setFixedPoints(vector<unsigned int>*);
	void setOutForce(PlainGeneralizedVector<Scalar>*);
	void getStaticPosition(vector<Vector<Scalar, Dim>> &);



private:
	vector<Vector<Scalar, Dim>> cur_pos_;
	ForceModelLinearSys<Scalar, Dim> linearSys_;
	FEMSolidForceModel<Scalar, Dim>* forcemodel_;
	VolumetricMesh<Scalar, Dim> *vMesh_;
	vector<unsigned int>* fixPoints_;
	PlainGeneralizedVector<Scalar>* force_;
	PlainGeneralizedVector<Scalar> dx_;
	PlainGeneralizedVector<Scalar> df;
};

template <typename Scalar, int Dim>
FEMelastostatic<Scalar,Dim>::FEMelastostatic(VolumetricMesh<Scalar, Dim> *vMesh, FEMSolidForceModel<Scalar, Dim> *forceModel, vector<unsigned int>*fixPoints, PlainGeneralizedVector<Scalar>*force):
vMesh_(vMesh),cur_pos_(),linearSys_(forceModel,fixPoints,&cur_pos_),fixPoints_(fixPoints),force_(force),forcemodel_(forceModel){
	dx_.resize(Dim*vMesh_->vertNum());
	df.resize(Dim*vMesh_->vertNum());
	for (unsigned int i = 0; i < vMesh_->vertNum(); ++i){
		cur_pos_.push_back(vMesh->vertPos(i));  //对cur_pos进行初始化
		for (int j = 0; j < Dim; ++j){
			dx_[i * Dim + j] = 0;
		}
	}
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
		forcemodel_->computeGlobalInternalForces(cur_pos_, cur_force);
		for (unsigned int i = 0; i < cur_force.size(); ++i){
			for (unsigned int kk = 0; kk < Dim; ++kk){
				df[i*Dim + kk] = force[i*Dim + kk] - cur_force[i][kk];
			}
		}
		linearSys_.filter(df);
		double dfsq = linearSys_.innerProduct(df, df);
		cout << "df*df:" << dfsq << endl;
		if (dfsq < 1e-6) break;
		for (unsigned int i = 0; i < dx_.size(); ++i)dx_[i] = 0;
		solver.solve(linearSys_, df, dx_);
		cout << "PCG Iterations:" << solver.iterationsUsed() << endl;
		//cout << dx_ << endl;
		for (int i = 0; i < cur_pos_.size(); ++i){
			for (unsigned int kk = 0; kk < Dim; ++kk){
				cur_pos_[i][kk] += dx_[i*Dim + kk];
			}
		}
		solver.reset();
		cout << "iteration :" << ++j << endl;
	} while (j < 100);

	staticPosition = cur_pos_;
}