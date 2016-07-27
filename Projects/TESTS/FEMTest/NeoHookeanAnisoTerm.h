/*
* @file  neo_hookeanAnisoTerm.h
* @brief Neo-Hookean hyperelastic material model with anisotropic term
* @author Fei Zhu
*
* This file is part of Physika, a versatile physics simulation library.
* Copyright (C) 2013- Physika Group.
*
* This Source Code Form is subject to the terms of the GNU General Public License v2.0.
* If a copy of the GPL was not distributed with this file, you can obtain one at:
* http://www.gnu.org/licenses/gpl-2.0.html
*
*/

#ifndef PHYSIKA_DYNAMICS_CONSTITUTIVE_MODELS_NEO_HOOKEAN_ANISO_TERM_H_
#define PHYSIKA_DYNAMICS_CONSTITUTIVE_MODELS_NEO_HOOKEAN_ANISO_TERM_H_

#include "Physika_Dynamics/Constitutive_Models/isotropic_hyperelastic_material.h"
#include "Physika_Dynamics/Constitutive_Models/constitutive_model_internal.h"
#include "Physika_Core/Vectors/vector_ND.h"

namespace Physika{

	template <typename Scalar, int Dim>
	class SquareMatrix;

	template <typename Scalar, int Dim>
	class NeoHookeanAnisoTerm : public IsotropicHyperelasticMaterial<Scalar, Dim>
	{
	public:
		NeoHookeanAnisoTerm();
		//if par_type = YOUNG_AND_POISSON, then: par1 = young's modulus, par2 = poisson_ratio
		//if par_type = LAME_COEFFICIENTS, then: par1 = lambda, par2 = mu
		NeoHookeanAnisoTerm(Scalar par1, Scalar par2, typename IsotropicHyperelasticMaterialInternal::ModulusType par_type,VectorND<Scalar> anisoDirection, Scalar aniso_C);
		NeoHookeanAnisoTerm(const NeoHookeanAnisoTerm<Scalar, Dim> &material);
		~NeoHookeanAnisoTerm();
		NeoHookeanAnisoTerm<Scalar, Dim>& operator= (const NeoHookeanAnisoTerm<Scalar, Dim> &material);
		NeoHookeanAnisoTerm<Scalar, Dim>* clone() const;
		void printInfo() const;
		Scalar energyDensity(const SquareMatrix<Scalar, Dim> &F) const;//compute potential energy density from given deformation gradient
		SquareMatrix<Scalar, Dim> firstPiolaKirchhoffStress(const SquareMatrix<Scalar, Dim> &F) const;
		SquareMatrix<Scalar, Dim> secondPiolaKirchhoffStress(const SquareMatrix<Scalar, Dim> &F) const;
		SquareMatrix<Scalar, Dim> cauchyStress(const SquareMatrix<Scalar, Dim> &F) const;
		//differential of first PiolaKirchhoff stress, for implicit time integration
		// \delta P = dP/dF : (\delta F)
		// \delta is differential
		virtual SquareMatrix<Scalar, Dim> firstPiolaKirchhoffStressDifferential(const SquareMatrix<Scalar, Dim> &F,
			const SquareMatrix<Scalar, Dim> &F_differential) const;

	protected:
		VectorND<Scalar> anisoDirection_;
		Scalar aniso_C_;
	};

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>::NeoHookeanAnisoTerm()
		:IsotropicHyperelasticMaterial<Scalar, Dim>()
	{
			anisoDirection_ = new VectorND<Scalar>(Dim, 0);
			aniso_C_ = 0;
	}

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>::NeoHookeanAnisoTerm(Scalar par1, Scalar par2, typename IsotropicHyperelasticMaterialInternal::ModulusType par_type, VectorND<Scalar> anisoDirection, Scalar aniso_C)
		: IsotropicHyperelasticMaterial<Scalar, Dim>(par1, par2, par_type)
	{
			anisoDirection_ = anisoDirection;
			aniso_C_ = aniso_C;
	}

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>::NeoHookeanAnisoTerm(const NeoHookeanAnisoTerm<Scalar, Dim> &material)
		: IsotropicHyperelasticMaterial<Scalar, Dim>(material)
	{
			anisoDirection_ = material.anisoDirection_;
			aniso_C_ = material.aniso_C_;
	}

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>::~NeoHookeanAnisoTerm()
	{
	}

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>& NeoHookeanAnisoTerm<Scalar, Dim>::operator= (const NeoHookeanAnisoTerm<Scalar, Dim> &material)
	{
		this->mu_ = material.mu_;
		this->lambda_ = material.lambda_;
		this->anisoDirection_ = material.anisoDirection_;
		this->aniso_C_ = material.aniso_C_;
		return *this;
	}

	template <typename Scalar, int Dim>
	NeoHookeanAnisoTerm<Scalar, Dim>* NeoHookeanAnisoTerm<Scalar, Dim>::clone() const
	{
		return new NeoHookeanAnisoTerm<Scalar, Dim>(*this);
	}

	template <typename Scalar, int Dim>
	void NeoHookeanAnisoTerm<Scalar, Dim>::printInfo() const
	{
		std::cout << "Compressible Neo-Hookean material with anisotropic term:" << std::endl;
		std::cout << "Energy density: Psi = mu/2*(trace(C)-Dim)-mu*lnJ+lambda/2*(lnJ)^2 +C*(|Fv|-1)^2" << std::endl;
		std::cout << "Where: C = transpose(F)*F, J = det(F)." << std::endl;
	}

	template <typename Scalar, int Dim>
	Scalar NeoHookeanAnisoTerm<Scalar, Dim>::energyDensity(const SquareMatrix<Scalar, Dim> &F) const
	{
		VectorND<Scalar> Fv(Dim, 0);
		for (unsigned int i = 0; i < Dim; ++i)
		for (unsigned int j = 0; j < Dim; ++j){
				Fv[j] += F(i, j)*anisoDirection_[j];
			}

		Scalar lenth = Fv.norm();
		Scalar trace_c = (F.transpose()*F).trace();
		Scalar J = F.determinant();
		Scalar lnJ = log(J);
		Scalar mu = this->mu_;
		Scalar lambda = this->lambda_;
		Scalar energy = mu / 2 * (trace_c - Dim) - mu*lnJ + lambda / 2 * lnJ*lnJ + aniso_C_*(lenth - 1)*(lenth - 1);
		return energy;
	}

	template <typename Scalar, int Dim>
	SquareMatrix<Scalar, Dim> NeoHookeanAnisoTerm<Scalar, Dim>::firstPiolaKirchhoffStress(const SquareMatrix<Scalar, Dim> &F) const
	{
		SquareMatrix<Scalar, Dim> identity = SquareMatrix<Scalar, Dim>::identityMatrix();
		SquareMatrix<Scalar, Dim> inverse_c = (F.transpose()*F).inverse();
		Scalar lnJ = log(F.determinant());
		Scalar mu = this->mu_;
		Scalar lambda = this->lambda_;
		SquareMatrix<Scalar, Dim> S = mu*(identity - inverse_c) + lambda*lnJ*inverse_c;
		VectorND<Scalar> Fv(Dim, 0);
		for (unsigned int i = 0; i < Dim; ++i)
		for (unsigned int j = 0; j < Dim; ++j){
			Fv[j] += F(i, j)*anisoDirection_[j];
		}
		Scalar lenth = Fv.norm();
		SquareMatrix<Scalar, Dim> FvvT;
		for (unsigned int i = 0; i < Dim;++i)
		for (unsigned int j = 0; j < Dim; ++j){
			FvvT(i, j) = Fv[i] * anisoDirection_[j];
		}

		SquareMatrix<Scalar, Dim> P = F*S + 2*aniso_C_*(lenth-1)/lenth*FvvT;
		return P;
	}

	template <typename Scalar, int Dim>
	SquareMatrix<Scalar, Dim> NeoHookeanAnisoTerm<Scalar, Dim>::secondPiolaKirchhoffStress(const SquareMatrix<Scalar, Dim> &F) const
	{
		SquareMatrix<Scalar, Dim> inverse_F = F.inverse();
		SquareMatrix<Scalar, Dim> S = inverse_F* firstPiolaKirchhoffStress(F);
		return S;
	}

	template <typename Scalar, int Dim>
	SquareMatrix<Scalar, Dim> NeoHookeanAnisoTerm<Scalar, Dim>::cauchyStress(const SquareMatrix<Scalar, Dim> &F) const
	{
		Scalar J = F.determinant();
		SquareMatrix<Scalar, Dim> stress = 1 / J*firstPiolaKirchhoffStress(F)*F.transpose();
		return stress;
	}

	template <typename Scalar, int Dim>
	SquareMatrix<Scalar, Dim> NeoHookeanAnisoTerm<Scalar, Dim>::firstPiolaKirchhoffStressDifferential(
		const SquareMatrix<Scalar, Dim> &F,
		const SquareMatrix<Scalar, Dim> &F_differential) const
	{
		VectorND<Scalar> Fv(Dim, 0);
		for (unsigned int i = 0; i < Dim; ++i)
		for (unsigned int j = 0; j < Dim; ++j){
			Fv[j] += F(i, j)*anisoDirection_[j];
		}
		Scalar lenth = Fv.norm();
		SquareMatrix<Scalar, Dim> FvvT;
		for (unsigned int i = 0; i < Dim; ++i)
		for (unsigned int j = 0; j < Dim; ++j){
			FvvT(i, j) = Fv[i] * anisoDirection_[j];
		}
		Scalar mu = this->mu_;
		Scalar lambda = this->lambda_;
		Scalar J = F.determinant();
		Scalar lnJ = log(J);
		SquareMatrix<Scalar, Dim> He = J*(F.inverse()).transpose();
		return mu*F_differential + ((lambda - (lambda*lnJ - mu)) / (J*J)*He.doubleContraction(F_differential))*He
			+ (lambda*lnJ - mu) / J*ConstitutiveModelInternal::cofactorMatrixDifferential(F, F_differential) + 2*aniso_C_*(F_differential*(1-1/lenth) + FvvT.doubleContraction(F_differential)/(lenth*lenth*lenth)*F);
	}

}  //end of namespace Physika

#endif //PHYSIKA_DYNAMICS_CONSTITUTIVE_MODELS_NEO_HOOKEAN_ANISO_TERM_H_
