#pragma once

#include "Materials/HyperelasticModel.h"
#include "Materials/MaterialPointSubstance.h"

namespace PhysX {

template <int Dim, Hyperelastic<Dim> Model = FixedCorotatedModel<Dim>>
class MaterialPointSoftBody : public MaterialPointSubstance<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using MaterialPointSubstance<Dim>::particles;
	using MaterialPointSubstance<Dim>::velocities;
	using MaterialPointSubstance<Dim>::velocityDerivatives;

protected:

	ParticlesBasedData<Dim, MatrixDr> _deformationGradients;

	const real _lameLambda;
	const real _lameMu;

public:

	MaterialPointSoftBody(const std::string &name, const Vector4f &color, const real density, const real lambda, const real mu) :
		MaterialPointSubstance<Dim>(name, color, density), _lameLambda(lambda), _lameMu(mu)
	{ }

	virtual ~MaterialPointSoftBody() = default;

	virtual void save(std::ofstream &fout) const override
	{
		MaterialPointSubstance<Dim>::save(fout);
		_deformationGradients.save(fout);
	}

	virtual void load(std::ifstream &fin) override
	{
		MaterialPointSubstance<Dim>::load(fin);
		_deformationGradients.load(fin);
	}

	virtual void reinitialize() override
	{
		MaterialPointSubstance<Dim>::reinitialize();
		_deformationGradients.resize(&particles);
		_deformationGradients.setConstant(MatrixDr::Identity());
	}

	virtual void update(const real dt) override
	{
		particles.parallelForEach([&](const int i) {
			particles.positions[i] += velocities[i] * dt;
			_deformationGradients[i] = (MatrixDr::Identity() + velocityDerivatives[i] * dt) * _deformationGradients[i];
		});
	}

	virtual void computeStressTensors(ParticlesBasedData<Dim, MatrixDr> &stresses) const override
	{
		stresses.resize(&particles);
		particles.parallelForEach([&](const int i) {
			stresses[i] = Model::computeStressTensorMultipliedByJ(_deformationGradients[i], _lameLambda, _lameMu);
		});
	}
};

}
