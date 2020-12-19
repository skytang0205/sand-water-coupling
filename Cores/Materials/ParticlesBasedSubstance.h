#pragma once

#include "Structures/ParticlesBasedData.h"

namespace PhysX {

template <int Dim>
class ParticlesBasedSubstance
{
	DECLARE_DIM_TYPES(Dim)

public:

	Particles<Dim> particles;
	ParticlesBasedVectorData<Dim> velocities;
	ParticlesBasedData<Dim, MatrixDr> velocityDerivatives;
	ParticlesBasedData<Dim, MatrixDr> deformationGradients;
	ParticlesBasedScalarData<Dim> plasticJacobians;

protected:

	const std::string _name;
	const Vector4f _color;

	const real _density;
	const real _lameLambda;
	const real _lameMu;

	const real _plasticLowerBound;
	const real _plasticUpperBound;
	const real _hardeningCoeff;

public:

	ParticlesBasedSubstance(const std::string &name, const Vector4f &color, const real density, const real lambda, const real mu, const real thetaC = -1, const real thetaS = -1, const real hardeningCoeff = 0) :
		_name(name),
		_color(color),
		_density(density),
		_lameLambda(lambda),
		_lameMu(mu),
		_plasticLowerBound(1 - thetaC),
		_plasticUpperBound(1 + thetaS),
		_hardeningCoeff(hardeningCoeff)
	{ }

	virtual ~ParticlesBasedSubstance() = default;

	const std::string &name() const { return _name; }
	const Vector4f &color() const { return _color; }
	real density() const { return _density; }
	bool plastic() const { return _plasticLowerBound <= _plasticUpperBound; }

	MatrixDr computeStressTensor(const int idx);
	void reinitialize();
};

}
