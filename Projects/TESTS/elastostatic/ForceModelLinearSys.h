#include "Physika_Dynamics/FEM/FEM_Solid_Force_Model/fem_solid_force_model.h"
#include "Physika_Core/Vectors/vector.h"
#include "Physika_Numerics/Linear_System_Solvers/linear_system.h"
#include "Physika_Numerics/Linear_System_Solvers/conjugate_gradient_solver.h"
#include "Physika_Numerics/Linear_System_Solvers/plain_generalized_vector_T.h"
#include<vector>

using std::vector;
using Physika::Vector;
using Physika::LinearSystem;
using Physika::GeneralizedVector;
using Physika::PlainGeneralizedVector;

template <typename Scalar, int Dim>
class ForceModelLinearSys : public LinearSystem<Scalar>{
public:
	ForceModelLinearSys(Physika::FEMSolidForceModel<Scalar, Dim> *forcemodel, vector<unsigned int> *fixedPoints, vector<Vector<double, 3>> *cur_pos){
		forcemodel_ = forcemodel;
		fixedPoints_ = fixedPoints;
		cur_pos_ = cur_pos;
	}
	virtual ~ForceModelLinearSys(){

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
	Physika::FEMSolidForceModel<Scalar, Dim> * forcemodel_;
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