#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

#include <concepts>

namespace PhysX {

template <typename Model, int Dim>
concept Hyperelastic = requires (const Matrix<Dim, real> &F, const real lambda, const real mu) {
	{ Model::computeNominalStressTensor(F, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeStressTensorMultipliedByJ(F, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
	{ Model::computeDeltaNominalStressTensor(F, F, lambda, mu) } -> std::convertible_to<Matrix<Dim, real>>;
};

template <int Dim>
class LinearModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static MatrixDr computeNominalStressTensor(const MatrixDr &F, const real lambda, const real mu)
	{
		return mu * (F + F.transpose() - 2 * MatrixDr::Identity()) + lambda * (F.trace() - Dim) * MatrixDr::Identity();
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &F, const real lambda, const real mu)
	{
		return computeNominalStressTensor(F, lambda, mu) * F.transpose();
	}

	static MatrixDr computeDeltaNominalStressTensor(const MatrixDr &F, const MatrixDr &dF, const real lambda, const real mu)
	{
		return mu * (dF + dF.transpose()) + lambda * dF.trace() * MatrixDr::Identity();
	}
};

template <int Dim>
class StVenantKirchhoffModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static MatrixDr computeNominalStressTensor(const MatrixDr &F, const real lambda, const real mu)
	{
		const MatrixDr FFT = F * F.transpose();
		return mu * (FFT - MatrixDr::Identity()) * F + lambda / 2 * (FFT.trace() - Dim) * F;
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &F, const real lambda, const real mu)
	{
		const MatrixDr FFT = F * F.transpose();
		return mu * (FFT - MatrixDr::Identity()) * FFT + lambda / 2 * (FFT.trace() - Dim) * FFT;
	}

	static MatrixDr computeDeltaNominalStressTensor(const MatrixDr &F, const MatrixDr &dF, const real lambda, const real mu)
	{
		const MatrixDr FTF = F.transpose() * F;
		const MatrixDr dFTF = dF.transpose() * F;
		return dF * mu * FTF + lambda / 2 * (FTF.trace() - Dim) * dF + F * mu * (dFTF + dFTF.transpose()) + lambda * dFTF.trace() * F;
	}
};

template <int Dim>
class CorotatedLinearModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static MatrixDr computeNominalStressTensor(const MatrixDr &F, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const MatrixDr R = svd.matrixU() * svd.matrixV().transpose();
		return 2 * mu * (F - R) + lambda * ((R * F.transpose()).trace() - Dim) * R;
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &F, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const MatrixDr RFT = svd.matrixU() * svd.matrixV().transpose() * F.transpose();
		return 2 * mu * (F * F.transpose() - RFT) + lambda * (RFT.trace() - Dim) * RFT;
	}

	static MatrixDr computeDeltaNominalStressTensor(const MatrixDr &F, const MatrixDr &dF, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const MatrixDr R = svd.matrixU() * svd.matrixV();
		const MatrixDr S = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
		// TODO: A dimension-free solution.
		const MatrixDr RTdFmdFTR = R.transpose() * dF - dF.transpose() * R;
		MatrixDr dR;
		if constexpr (Dim == 2) {
			const real x = RTdFmdFTR(0, 1) / S.trace();
			dR << 0, x, -x, 0;
		}
		else {
			MatrixDr A;
			A << S(0, 0) + S(1, 1), S(1, 2), -S(0, 2), S(1, 2), S(0, 0) + S(2, 2), S(0, 1), -S(0, 2), S(0, 1), S(1, 1) + S(2, 2);
			VectorDr x = A.inverse() * VectorDr(RTdFmdFTR(0, 1), RTdFmdFTR(0, 2), RTdFmdFTR(1, 2));
			dR << 0, x[0], x[1], -x[0], 0, x[2], -x[1], -x[2], 0;
		}
		dR = R * dR;
		return 2 * mu * (dF - dR) + lambda * ((R * F.transpose()).trace() - Dim) * dR + lambda * (dR * F.transpose() + R * dF.transpose()).trace() * R;
	}
};


template <int Dim>
class FixedCorotatedLinearModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static MatrixDr computeNominalStressTensor(const MatrixDr &F, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real J = svd.singularValues().prod();
		return 2 * mu * (F - svd.matrixU() * svd.matrixV().transpose()) + lambda * J * (J - 1) * F.transpose().inverse();
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &F, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real J = svd.singularValues().prod();
		return 2 * mu * (F - svd.matrixU() * svd.matrixV().transpose()) * F.transpose() + lambda * J * (J - 1) * MatrixDr::Identity();
	}

	static MatrixDr computeDeltaNominalStressTensor(const MatrixDr &F, const MatrixDr &dF, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const real J = svd.singularValues().prod();
		const MatrixDr R = svd.matrixU() * svd.matrixV();
		const MatrixDr S = svd.matrixV() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
		const MatrixDr invF = F.inverse();
		// TODO: A dimension-free solution.
		const MatrixDr RTdFmdFTR = R.transpose() * dF - dF.transpose() * R;
		MatrixDr dR;
		if constexpr (Dim == 2) {
			const real x = RTdFmdFTR(0, 1) / S.trace();
			dR << 0, x, -x, 0;
		}
		else {
			MatrixDr A;
			A << S(0, 0) + S(1, 1), S(1, 2), -S(0, 2), S(1, 2), S(0, 0) + S(2, 2), S(0, 1), -S(0, 2), S(0, 1), S(1, 1) + S(2, 2);
			VectorDr x = A.inverse() * VectorDr(RTdFmdFTR(0, 1), RTdFmdFTR(0, 2), RTdFmdFTR(1, 2));
			dR << 0, x[0], x[1], -x[0], 0, x[2], -x[1], -x[2], 0;
		}
		dR = R * dR;
		return 2 * mu * (dF - dR) + lambda * J * invF.transpose() * (J * invF.transpose().cwiseProduct(dF)) + lambda * J * (J - 1) * ((invF * dF).trace() * invF - invF * dF * invF).transpose();
	}
};

template <int Dim>
class NeoHookeanModel
{
	DECLARE_DIM_TYPES(Dim)

public:

	static MatrixDr computeNominalStressTensor(const MatrixDr &F, const real lambda, const real mu)
	{
		const real J = F.determinant();
		const MatrixDr invF = F.inverse();
		return mu * (F - invF.transpose()) + lambda * std::log(J) * invF.transpose();
	}

	static MatrixDr computeStressTensorMultipliedByJ(const MatrixDr &F, const real lambda, const real mu)
	{
		const real J = F.determinant();
		return mu * (F * F.transpose() - MatrixDr::Identity()) + lambda * std::log(J) * MatrixDr::Identity();
	}

	static MatrixDr computeDeltaNominalStressTensor(const MatrixDr &F, const MatrixDr &dF, const real lambda, const real mu)
	{
		const real J = F.determinant();
		const MatrixDr invF = F.inverse();
		return mu * dF + ((mu - lambda * std::log(J)) * invF * dF * invF + lambda * (invF * dF).trace() * invF).transpose();
	}
};

}
