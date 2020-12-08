#pragma once

#include "Structures/Particles.h"
#include "Structures/ParticlesNearbySearcher.h"

namespace PhysX {

template <int Dim>
class SmoothedParticles final : public Particles<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using Particles<Dim>::positions;

protected:

	using Particles<Dim>::_mass;

	const real _radius;
	const real _squaredRadius;
	const real _invRadius;
	const real _invSquaredRadius;

	const real _stdKernelNormCoeff0;
	const real _stdKernelNormCoeff1;

	const real _spikyKernelNormCoeff0;
	const real _spikyKernelNormCoeff1;

	std::unique_ptr<ParticlesNearbySearcher<Dim>> _nearbySearcher;

public:

	SmoothedParticles(const real mass, const real radius, const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero());

	SmoothedParticles &operator=(const SmoothedParticles &rhs) = delete;
	virtual ~SmoothedParticles() = default;

	real radius() const { return _radius; }

	// Interpolation helper functions.

	real stdKernel(const real distance) const;
	real firstDerivativeStdKernel(const real distance) const;
	VectorDr gradientStdKernel(const real distance, const VectorDr &direction) const;
	real secondDerivativeStdKernel(const real distance) const;

	real spikyKernel(const real distance) const;
	real firstDerivativeSpikyKernel(const real distance) const;
	VectorDr gradientSpikyKernel(const real distance, const VectorDr &direction) const;
	real secondDerivativeSpikyKernel(const real distance) const;

	void resetNeighborSearcher() { _nearbySearcher->reset(positions); }
	void forEachNearby(const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func) const { _nearbySearcher->forEach(positions, pos, func); }

protected:
};

}
