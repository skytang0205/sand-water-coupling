#pragma once

#include "MaterialPointSoftBody.h"

namespace PhysX {

template <int Dim, Hyperelastic<Dim> Model = FixedCorotatedModel<Dim>>
class MatPointPlasticSoftBody : public MaterialPointSoftBody<Dim, Model>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using MaterialPointSubstance<Dim>::particles;
	using MaterialPointSubstance<Dim>::velocities;
	using MaterialPointSubstance<Dim>::velocityDerivatives;

protected:

	using MaterialPointSoftBody<Dim, Model>::_deformationGradients;
	using MaterialPointSoftBody<Dim, Model>::_lameLambda;
	using MaterialPointSoftBody<Dim, Model>::_lameMu;

	ParticlesBasedScalarData<Dim> _plasticJacobians;

	const real _plasticLowerBound;
	const real _plasticUpperBound;
	const real _hardeningCoeff;

public:

	MatPointPlasticSoftBody(const std::string &name, const Vector4f &color, const real density, const real lambda, const real mu, const real thetaC, const real thetaS, const real hardeningCoeff) :
		MaterialPointSoftBody<Dim, Model>(name, color, density, lambda, mu),
		_plasticLowerBound(1 - thetaC),
		_plasticUpperBound(1 + thetaS),
		_hardeningCoeff(hardeningCoeff)
	{ }

	virtual ~MatPointPlasticSoftBody() = default;

	virtual void save(std::ofstream &fout) const override
	{
		MaterialPointSoftBody<Dim, Model>::save(fout);
		_plasticJacobians.save(fout);
	}

	virtual void load(std::ifstream &fin) override
	{
		MaterialPointSoftBody<Dim, Model>::load(fin);
		_plasticJacobians.load(fin);
	}

	virtual void reinitialize() override
	{
		MaterialPointSoftBody<Dim, Model>::reinitialize();
		_plasticJacobians.resize(&particles);
		_plasticJacobians.setConstant(1);
	}

	virtual void update(const real dt) override
	{
		particles.parallelForEach([&](const int i) {
			particles.positions[i] += velocities[i] * dt;
			_deformationGradients[i] = (MatrixDr::Identity() + velocityDerivatives[i] * dt) * _deformationGradients[i];

			Eigen::JacobiSVD<MatrixDr> svd(_deformationGradients[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
			_plasticJacobians[i] *= svd.singularValues().prod();
			const VectorDr svdS = svd.singularValues().cwiseMax(_plasticLowerBound).cwiseMin(_plasticUpperBound);
			_plasticJacobians[i] /= svdS.prod();
			_deformationGradients[i] = svd.matrixU() * svdS.asDiagonal() * svd.matrixV().transpose();
		});
	}

	virtual void computeStressTensors(ParticlesBasedData<Dim, MatrixDr> &stresses) const override
	{
		stresses.resize(&particles);
		particles.parallelForEach([&](const int i) {
			const real ratio = std::exp(_hardeningCoeff * (1 - _plasticJacobians[i]));
			stresses[i] = Model::computeStressTensorMultipliedByJ(_deformationGradients[i], ratio * _lameLambda, ratio * _lameMu);
		});
	}
};

}
