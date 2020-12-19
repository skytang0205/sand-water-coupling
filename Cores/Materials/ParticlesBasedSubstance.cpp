#include "ParticlesBasedSubstance.h"

namespace PhysX {

template <int Dim>
Matrix<Dim, real> ParticlesBasedSubstance<Dim>::computeStressTensor(int idx)
{
	real lambda = _lameLambda;
	real mu = _lameMu;
	if (plastic()) {
		const real ratio = std::exp(_hardeningCoeff * (1 - plasticJacobians[idx]));
		lambda *= ratio;
		mu *= ratio;
	}

	MatrixDr &matF = deformationGradients[idx];

	Eigen::JacobiSVD<MatrixDr> svd(matF, Eigen::ComputeFullU | Eigen::ComputeFullV);
	VectorDr svdS = svd.singularValues();
	const MatrixDr svdU = svd.matrixU();
	const MatrixDr svdV = svd.matrixV();

	if (plastic()) {
		// Handle plasticity.
		plasticJacobians[idx] *= svdS.prod();
		svdS = svdS.cwiseMax(_plasticLowerBound).cwiseMin(_plasticUpperBound);
	}
	const real jacobian = svdS.prod();

	if (plastic()) {
		plasticJacobians[idx] /= jacobian;
		matF = svdU * svdS.asDiagonal() * svdV.transpose();
	}

	if (!mu) {
		// For fluid, reset deformation gradient to avoid numerical instability.
		if constexpr (Dim == 2) matF = MatrixDr::Identity() * std::sqrt(jacobian);
		else matF = MatrixDr::Identity() * std::cbrt(jacobian);
	}

	// Compute by the fixed corotated constitutive model.
	// The result is actually P(F) * F^T
	return 2 * mu * (matF - svdU * svdV.transpose()) * matF.transpose() + MatrixDr::Identity() * lambda * jacobian * (jacobian - 1);
}

template <int Dim>
void ParticlesBasedSubstance<Dim>::reinitialize()
{
	velocities.resize(&particles);
	velocities.setZero();

	velocityDerivatives.resize(&particles);
	velocityDerivatives.setZero();

	deformationGradients.resize(&particles);
	deformationGradients.setConstant(MatrixDr::Identity());

	if (plastic()) {
		plasticJacobians.resize(&particles);
		plasticJacobians.setConstant(1);
	}
}

template class ParticlesBasedSubstance<2>;
template class ParticlesBasedSubstance<3>;

}
