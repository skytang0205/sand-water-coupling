#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

namespace PhysX {

template <typename Model, int Dim>
concept Hyperelastic = requires (const Matrix<Dim, real> &defmGrad, const real lambda, const real mu) {
	Model::energy(defmGrad, lambda, mu);
	Model::nominalStressTensor(defmGrad, lambda, mu);
};

template <int Dim>
class NeoHookeanModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static real energy(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real logJ = std::log(defmGrad.determinant());
		return mu / 2 * ((defmGrad * defmGrad.transpose()).trace() - Dim) - mu * logJ + lambda / 2 * logJ * logJ;
	}

	static MatrixDr nominalStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		const real logJ = std::log(defmGrad.determinant());
		const MatrixDr invDefmGradT = defmGrad.transpose().inverse();
		return mu * (defmGrad - invDefmGradT) + lambda * logJ * invDefmGradT;
	}

	static MatrixDr trueStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
	}
};

template <int Dim>
class FixedCorotatedModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static real energy(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad);
		return (svd.singularValues() - VectorDr::Ones()).squaredNorm() * mu + MathFunc::square(svd.singularValues().prod() - 1) * lambda / 2;
	}

	static MatrixDr nominalStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real jacobian = svd.singularValues().prod();
		2 * mu * (defmGrad - svd.matrixU() * svd.matrixV().transpose()) + lambda * jacobian * (jacobian - 1) * defmGrad.transpose().inverse();
	}

	static MatrixDr trueStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
	}
};

}
