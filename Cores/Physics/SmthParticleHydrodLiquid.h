#pragma once

#include "Physics/Simulation.h"
#include "Geometries/Collider.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"

namespace PhysX {

template <int Dim>
class SmthParticleHydrodLiquid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class SmthPartHydrodLiquidBuilder;

protected:

	const real _particleSpacing;

	SmoothedParticles<Dim> _particles;
	ParticlesBasedVectorField<Dim> _velocities;
	ParticlesBasedScalarField<Dim> _pressures;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

	bool _enableGravity = true;

	real _viscosityCoeff = real(.01);
	real _pseudoViscosityCoeff = 10;
	real _targetDensity = real(1e3);
	real _speedOfSound = 1482;
	real _eosExponent = 7;

public:

	SmthParticleHydrodLiquid(const real particleMass, const real particleSpacing) :
		_particleSpacing(particleSpacing),
		_particles(particleMass, _particleSpacing * 2)
	{ }

	SmthParticleHydrodLiquid(const SmthParticleHydrodLiquid &rhs) = delete;
	SmthParticleHydrodLiquid &operator=(const SmthParticleHydrodLiquid &rhs) = delete;
	virtual ~SmthParticleHydrodLiquid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return std::min(real(1) / frameRate, stepRate * _particles.radius() / _velocities.normMax()); }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	void reinitializeParticlesBasedData();
	void moveParticles(const real dt);
	void applyExternalForces(const real dt);
	void applyViscosityForce(const real dt);
	void applyPressureForce(const real dt);
	void applyPseudoViscosity(const real dt);
};

}
