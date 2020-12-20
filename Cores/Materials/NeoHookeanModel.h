#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

namespace PhysX {

template <int Dim>
class NeoHookeanModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static real computeEnergy(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real logJ = std::log(defmGrad.determinant());
		return mu / 2 * ((defmGrad * defmGrad.transpose()).trace() - Dim) - mu * logJ + lambda / 2 * logJ * logJ;
	}

	static MatrixDr computeNominalStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real jacobian = defmGrad.determinant();
		const MatrixDr invDefmGradT = defmGrad.transpose().inverse();
		return mu * (defmGrad - invDefmGradT) + lambda * std::log(jacobian) * invDefmGradT;
	}

	static MatrixDr computeStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real jacobian = defmGrad.determinant();
		return (mu * (defmGrad * defmGrad.transpose() - MatrixDr::Identity()) + lambda * std::log(jacobian) * MatrixDr::Identity()) / jacobian;
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real jacobian = defmGrad.determinant();
		return mu * (defmGrad * defmGrad.transpose() - MatrixDr::Identity()) + lambda * std::log(jacobian) * MatrixDr::Identity();
	}

	static Matrix<Dim * Dim, real> computeEnergyHessian(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		// TODO: Impl
		return Matrix<Dim * Dim, real>::Zero();
	}
};

}
