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

	using MaterialPointSoftBody<Dim>::_deformationGradients;

	ParticlesBasedScalarData<Dim> _plasticJacobians;

	const real _plasticUpperBound;
	const real _plasticUpperBound;
	const real _hardeningCoeff;

public:

	MatPointPlasticSoftBody(const std::string &name, const Vector4f &color, const real density, const real lambda, const real mu, const real thetaC, const real thetaS, const real hardeningCoeff) :
		MaterialPointSoftBody<Dim>(name, color, density, _lameLambda, _lameMu)
		_plasticLowerBound(1 - thetaC),
		_plasticUpperBound(1 + thetaS),
		_hardeningCoeff(hardeningCoeff)
	{ }

	virtual ~MatPointPlasticSoftBody() = default;
};

}
