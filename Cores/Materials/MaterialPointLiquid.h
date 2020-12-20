#pragma once

#include "MaterialPointSubstance.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class MaterialPointLiquid : public MaterialPointSubstance<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using MaterialPointSubstance<Dim>::particles;
	using MaterialPointSubstance<Dim>::velocities;
	using MaterialPointSubstance<Dim>::velocityDerivatives;

protected:

	ParticlesBasedScalarData<Dim> _jacobians;

	const real _bulkModulus;

public:

	MaterialPointLiquid(const std::string &name, const Vector4f &color, const real density, const real bulkModulus) :
		MaterialPointSubstance<Dim>(name, color, density),
		_bulkModulus(bulkModulus)
	{ }

	virtual ~MaterialPointLiquid() = default;

	virtual void save(std::ofstream &fout) const override
	{
		MaterialPointSubstance<Dim>::save(fout);
		_jacobians.save(fout);
	}

	virtual void load(std::ifstream &fin) override
	{
		MaterialPointSubstance<Dim>::load(fin);
		_jacobians.load(fin);
	}

	virtual void reinitialize() override
	{
		MaterialPointSubstance<Dim>::reinitialize();
		_jacobians.resize(&particles);
		_jacobians.setConstant(1);
	}

	virtual void update(const real dt) override
	{
		particles.parallelForEach([&](const int i) {
			particles.positions[i] += velocities[i] * dt;
			_jacobians[i] *= (1 + velocityDerivatives[i].trace() * dt);
		});
	}

	virtual ParticlesBasedData<Dim, MatrixDr> getDeformationGradients() const override
	{
		ParticlesBasedData<Dim, MatrixDr> defmGrads(&particles);
		particles.parallelForEach([&](const int i) {
			if constexpr (Dim == 2)
				defmGrads[i] = MatrixDr::Identity() * std::sqrt(_jacobians[i]);
			else
				defmGrads[i] = MatrixDr::Identity() * std::cbrt(_jacobians[i]);
		});
		return defmGrads;
	}

	virtual void computeStressTensors(ParticlesBasedData<Dim, MatrixDr> &stresses) const override
	{
		stresses.resize(&particles);
		particles.parallelForEach([&](const int i) {
			stresses[i] = MatrixDr::Identity() * (_jacobians[i] - 1) * _bulkModulus * _jacobians[i];
		});
	}

	virtual void computeEnergyHessians(ParticlesBasedData<Dim, Matrix<Dim * Dim, real>> &hessians) const override
	{
		hessians.resize(&particles);
		particles.parallelForEach([&](const int i) {
			// TODO: Impl
		});
	}
};

}
