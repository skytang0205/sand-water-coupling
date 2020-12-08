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

	ParticlesScalarAttribute<Dim> densities;

protected:

	using Particles<Dim>::_mass;

	const real _radius;
	const real _squaredRadius;
	const real _invRadius;
	const real _invSquaredRadius;

	const real _stdKernelNormCoeff0;
	const real _stdKernelNormCoeff1;
	const real _stdKernelNormCoeff2;

	const real _spikyKernelNormCoeff0;
	const real _spikyKernelNormCoeff1;
	const real _spikyKernelNormCoeff2;

	std::unique_ptr<ParticlesNearbySearcher<Dim>> _nearbySearcher;

public:

	SmoothedParticles(const real mass, const real radius, const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero());

	SmoothedParticles &operator=(const SmoothedParticles &rhs) = delete;
	virtual ~SmoothedParticles() = default;

	real radius() const { return _radius; }

	// Interpolation helper functions.

	real stdKernel(const real distance) const
	{
		if (const real x = 1 - distance * distance * _invSquaredRadius; x > 0)
			return _stdKernelNormCoeff0 * x * x * x;
		else return 0;
	}

	real laplacianStdKernel(const real distance) const
	{
		if (const real x = 1 - distance * distance * _invSquaredRadius; x > 0)
			return _stdKernelNormCoeff2 * x * ((Dim + 4) * x - 4);
		else return 0;
	}

	real stdKernel(const VectorDr &deltaPos) const { return stdKernel(deltaPos.norm()); }
	VectorDr gradientStdKernel(const VectorDr &deltaPos) const { return firstDerivativeStdKernel(deltaPos.norm()) * deltaPos.normalized(); }
	real laplacianStdKernel(const VectorDr &deltaPos) const { return laplacianStdKernel(deltaPos.norm()); }

	real spikyKernel(const real distance) const
	{
		if (const real x = 1 - distance * _invRadius; x > 0)
			return _spikyKernelNormCoeff0 * x * x * x;
		else return 0;
	}

	real spikyKernel(const VectorDr &deltaPos) const { return spikyKernel(deltaPos.norm()); }
	VectorDr gradientSpikyKernel(const VectorDr &deltaPos) const { return firstDerivativeSpikyKernel(deltaPos.norm()) * deltaPos.normalized(); }

	real laplacianSpikyKernel(const VectorDr &deltaPos) const
	{
		const real distance = deltaPos.norm();
		return secondDerivativeSpikyKernel(distance) + firstDerivativeSpikyKernel(distance) * (Dim - 1) / distance;
	}

	void resetNeighborSearcher() { _nearbySearcher->reset(positions); }
	void forEachNearby(const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func) const { _nearbySearcher->forEach(positions, pos, func); }

protected:

	real firstDerivativeStdKernel(const real distance) const
	{
		if (const real x = 1 - distance * distance * _invSquaredRadius; x > 0)
			return _stdKernelNormCoeff1 * x * x * distance;
		else return 0;
	}

	real secondDerivativeStdKernel(const real distance) const
	{
		if (const real x = 1 - distance * distance * _invSquaredRadius; x > 0)
			return _stdKernelNormCoeff2 * x * (5 * x - 4);
		else return 0;
	}

	real firstDerivativeSpikyKernel(const real distance) const
	{
		if (const real x = 1 - distance * _invRadius; x > 0)
			return _spikyKernelNormCoeff1 * x * x;
		else return 0;
	}

	real secondDerivativeSpikyKernel(const real distance) const
	{
		if (const real x = 1 - distance * _invRadius; x > 0)
			return _spikyKernelNormCoeff2 * x;
		else return 0;
	}
};

}
