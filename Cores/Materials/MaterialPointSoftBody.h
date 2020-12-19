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
	{
	}

	virtual ~MaterialPointSoftBody() = default;

	virtual void save(std::ofstream &fout) const;
	virtual void load(std::ifstream &fin);

	virtual void reinitialize() override;
	virtual void update(const real dt);

protected:

	bool plastic() const { return _plasticLowerBound <= _plasticUpperBound; }
};

}
