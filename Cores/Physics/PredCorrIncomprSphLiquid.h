#pragma once

#include "Physics/SmthParticleHydrodLiquid.h"

namespace PhysX {

template <int Dim>
class PredCorrIncomprSphLiquid : public SmthParticleHydrodLiquid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using SmthParticleHydrodLiquid<Dim>::_particles;
	using SmthParticleHydrodLiquid<Dim>::_velocities;
	using SmthParticleHydrodLiquid<Dim>::_pressures;
	using SmthParticleHydrodLiquid<Dim>::_colliders;
	using SmthParticleHydrodLiquid<Dim>::_targetDensity;

	static constexpr int _kPredCorrMaxIters = 5;
	static constexpr real _kPredCorrErrorRatio = real(.01);

	ParticlesBasedVectorData<Dim> _predPositions;
	ParticlesBasedVectorField<Dim> _predVelocities;
	ParticlesBasedScalarData<Dim> _densityErrors;

public:

	PredCorrIncomprSphLiquid(const real particleRadius) : SmthParticleHydrodLiquid<Dim>(particleRadius) { }

	PredCorrIncomprSphLiquid(const PredCorrIncomprSphLiquid &rhs) = delete;
	PredCorrIncomprSphLiquid &operator=(const PredCorrIncomprSphLiquid &rhs) = delete;
	virtual ~PredCorrIncomprSphLiquid() = default;

	virtual void advance(const real dt) override;

protected:

	using SmthParticleHydrodLiquid<Dim>::moveParticles;
	using SmthParticleHydrodLiquid<Dim>::applyExternalForces;
	using SmthParticleHydrodLiquid<Dim>::applyViscosityForce;

	virtual void reinitializeParticlesBasedData() override;
	virtual void applyPressureForce(const real dt) override;

	real computeDelta(const real dt) const;
};

}
