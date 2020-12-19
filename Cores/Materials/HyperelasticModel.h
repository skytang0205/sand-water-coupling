#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

#include <concepts>

namespace PhysX {

template <typename Model, int Dim>
concept Hyperelastic = requires (const Matrix<Dim, real> &defmGrad, const real lambda, const real mu) {
	{ Model::computeEnergy(defmGrad, lambda, mu) } -> std::convertible_to<real>;
	{ Model::computeNominalStressTensor(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeStressTensor(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeStressTensorMultipliedByJ(defmGrad, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
};

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
};

template <int Dim>
class FixedCorotatedModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static real computeEnergy(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad);
		return (svd.singularValues() - VectorDr::Ones()).squaredNorm() * mu + MathFunc::square(svd.singularValues().prod() - 1) * lambda / 2;
	}

	static MatrixDr computeNominalStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real jacobian = svd.singularValues().prod();
		return 2 * mu * (defmGrad - svd.matrixU() * svd.matrixV().transpose()) + lambda * jacobian * (jacobian - 1) * defmGrad.transpose().inverse();
	}

	static MatrixDr computeStressTensor(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real jacobian = svd.singularValues().prod();
		return (2 * mu * (defmGrad - svd.matrixU() * svd.matrixV().transpose()) * defmGrad.transpose() + lambda * jacobian * (jacobian - 1) * MatrixDr::Identity()) / jacobian;
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real jacobian = svd.singularValues().prod();
		return 2 * mu * (defmGrad - svd.matrixU() * svd.matrixV().transpose()) * defmGrad.transpose() + lambda * jacobian * (jacobian - 1) * MatrixDr::Identity();
	}
};

}
