#pragma once

#include "Physics/Simulation.h"
#include "Structures/ParticlesBasedScalarField.h"
#include "Structures/ParticlesBasedVectorField.h"

namespace PhysX {

template <int Dim>
class SmthParticleHydrodLiquid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

	SmoothedParticles<Dim> _particles;
	ParticlesBasedVectorField<Dim> _velocities;
	ParticlesBasedScalarField<Dim> _pressures;

	std::vector<std::unique_ptr<Collider<Dim>>> _colliders;

	bool _enableGravity = true;

	real _targetDensity = real(1e3);
	real _speedOfSound = 1482;
	real _eosExponent = 7;

public:

	SmthParticleHydrodLiquid() { }

	SmthParticleHydrodLiquid(const SmthParticleHydrodLiquid &rhs) = delete;
	SmthParticleHydrodLiquid &operator=(const SmthParticleHydrodLiquid &rhs) = delete;
	virtual ~SmthParticleHydrodLiquid() = default;

	virtual real getTimeStep(const uint frameRate, const real stepRate) const override { return 1; }

	virtual int dimension() const override { return Dim; }
	virtual void writeDescription(YAML::Node &root) const override;
	virtual void writeFrame(const std::string &frameDir, const bool staticDraw) const override;
	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

	virtual void initialize() override;
	virtual void advance(const real dt) override;

protected:

	void moveParticles(const real dt);
	void applyPressureForce(const real dt);
	void applyExternalForces(const real dt);
};

}
