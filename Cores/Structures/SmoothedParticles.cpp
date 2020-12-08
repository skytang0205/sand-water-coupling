#include "SmoothedParticles.h"

#include <numbers>

namespace PhysX {

template <int Dim>
SmoothedParticles<Dim>::SmoothedParticles(const real mass, const real radius, const size_t cnt, const VectorDr &pos) :
	Particles<Dim>(mass, cnt, pos),
	_radius(radius),
	_squaredRadius(_radius *_radius),
	_invRadius(1 / _radius),
	_invSquaredRadius(1 / _squaredRadius),
	_stdKernelNormCoeff0((Dim == 2 ? real(4) : 315 * _invRadius / 64) * real(std::numbers::inv_pi) * _invSquaredRadius),
	_stdKernelNormCoeff1(-6 * _invSquaredRadius * _stdKernelNormCoeff0),
	_stdKernelNormCoeff2(_stdKernelNormCoeff1),
	_spikyKernelNormCoeff0((Dim == 2 ? real(10) : 15 * _invRadius) * real(std::numbers::inv_pi) *_invSquaredRadius),
	_spikyKernelNormCoeff1(-3 * _invRadius * _spikyKernelNormCoeff0),
	_spikyKernelNormCoeff2(-2 * _invRadius * _spikyKernelNormCoeff1),
	_nearbySearcher(std::make_unique<ParticlesNearbySearcher<Dim>>(_radius))
{ }

template class SmoothedParticles<2>;
template class SmoothedParticles<3>;

}
