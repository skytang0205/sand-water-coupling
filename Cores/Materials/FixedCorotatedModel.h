#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

namespace PhysX {

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

	static Matrix<Dim * Dim, real> computeEnergyHessian(const MatrixDr &defmGrad, const real lambda, const real mu)
	{
		Eigen::JacobiSVD<MatrixDr> svd(defmGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
		const VectorDr svdS = svd.singularValues();
		const MatrixDr svdU = svd.matrixU();
		const MatrixDr svdV = svd.matrixV();

		Matrix<Dim * Dim, real> hessianS = Matrix<Dim * Dim, real>::Zero();

		if constexpr (Dim == 2) {
			hessianS(0, 0) = lambda * svdS[1] * svdS[1] + mu * 2;
			hessianS(3, 3) = lambda * svdS[0] * svdS[0] + mu * 2;
			hessianS(1, 1) = mu * 2 * (1 - 1 / svdS.sum());
			hessianS(2, 2) = mu * 2 * (1 - 1 / svdS.sum());

			hessianS(0, 3) = hessianS(3, 0) = lambda * (2 * svdS.prod() - 1);
			hessianS(1, 2) = hessianS(2, 1) = lambda * (1 - svdS.prod()) + 2 * mu / svdS.sum();
		}
		else {
			hessianS(0, 0) = lambda * svdS[1] * svdS[1] * svdS[2] * svdS[2] + mu * 2;
			hessianS(4, 4) = lambda * svdS[0] * svdS[0] * svdS[2] * svdS[2] + mu * 2;
			hessianS(8, 8) = lambda * svdS[0] * svdS[0] * svdS[1] * svdS[1] + mu * 2;
			hessianS(1, 1) = mu * 2 * (1 - 1 / (svdS[0] + svdS[1]));
			hessianS(3, 3) = mu * 2 * (1 - 1 / (svdS[0] + svdS[1]));
			hessianS(2, 2) = mu * 2 * (1 - 1 / (svdS[0] + svdS[2]));
			hessianS(6, 6) = mu * 2 * (1 - 1 / (svdS[0] + svdS[2]));
			hessianS(5, 5) = mu * 2 * (1 - 1 / (svdS[1] + svdS[2]));
			hessianS(7, 7) = mu * 2 * (1 - 1 / (svdS[1] + svdS[2]));

			hessianS(0, 4) = hessianS(4, 0) = lambda * svdS[2] * (2 * svdS.prod() - 1);
			hessianS(0, 8) = hessianS(8, 0) = lambda * svdS[1] * (2 * svdS.prod() - 1);
			hessianS(4, 8) = hessianS(8, 4) = lambda * svdS[0] * (2 * svdS.prod() - 1);
			hessianS(1, 3) = hessianS(3, 1) = lambda * svdS[2] * (1 - svdS.prod()) + 2 * mu / (svdS[0] + svdS[1]);
			hessianS(2, 6) = hessianS(6, 2) = lambda * svdS[1] * (1 - svdS.prod()) + 2 * mu / (svdS[0] + svdS[2]);
			hessianS(5, 7) = hessianS(7, 5) = lambda * svdS[0] * (1 - svdS.prod()) + 2 * mu / (svdS[1] + svdS[2]);
		}

		Matrix<Dim * Dim, real> hessian;
		for (int i = 0; i < Dim; i++)
			for (int j = 0; j < Dim; j++)
				for (int r = 0; r < Dim; r++)
					for (int s = 0; s < Dim; s++) {
						auto &element = hessian(i * Dim + j, r * Dim + s);
						element = 0;
						for (int k = 0; k < Dim; k++)
							for (int l = 0; l < Dim; l++)
								for (int m = 0; m < Dim; m++)
									for (int n = 0; n < Dim; n++)
										element += hessianS(k * Dim + l, m * Dim + n) * svdU(i, k) * svdU(r, m) * svdV(s, n) * svdV(j, l);
					}
		return hessian;
	}
};

}
