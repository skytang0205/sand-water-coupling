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

	virtual void update(const int idx, const real dt) override
	{
		MaterialPointSubstance<Dim>::update(idx, dt);
		_jacobians[idx] *= (1 + velocityDerivatives[idx].trace() * dt);
	}

	virtual MatrixDr computeStressTensor(const int idx) const override
	{
		return MatrixDr::Identity() * (_jacobians[idx] - 1) * _bulkModulus * _jacobians[idx];
	}

	virtual MatrixDr computeDeltaStressTensor(const int idx, const MatrixDr &weightSum) const override
	{
		return MatrixDr::Zero();
	}
};

}
