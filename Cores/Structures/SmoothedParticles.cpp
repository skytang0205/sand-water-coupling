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
	_spikyKernelNormCoeff0((Dim == 2 ? real(10) : 15 * _invRadius) * real(std::numbers::inv_pi) *_invSquaredRadius),
	_spikyKernelNormCoeff1(-3 * _invRadius * _spikyKernelNormCoeff0),
	_nearbySearcher(std::make_unique<ParticlesNearbySearcher<Dim>>(_radius))
{ }

template <int Dim>
real SmoothedParticles<Dim>::stdKernel(const real distance) const
{
	if (distance * distance < _squaredRadius) {
		const real x = 1 - distance * distance * _invSquaredRadius;
		return _stdKernelNormCoeff0 * x * x * x;
	}
	else return 0;
}

template <int Dim>
real SmoothedParticles<Dim>::firstDerivativeStdKernel(const real distance) const
{
	if (distance * distance < _squaredRadius) {
		const real x = 1 - distance * distance * _invSquaredRadius;
		return _stdKernelNormCoeff1 * distance * x * x;
	}
	else return 0;
}

template <int Dim>
real SmoothedParticles<Dim>::spikyKernel(const real distance) const
{
	if (distance * distance < _squaredRadius) {
		const real x = 1 - distance * _invRadius;
		return _spikyKernelNormCoeff0 * x * x * x;
	}
	else return 0;
}

template <int Dim>
real SmoothedParticles<Dim>::firstDerivativeSpikyKernel(const real distance) const
{
	if (distance * distance < _squaredRadius) {
		const real x = 1 - distance * _invRadius;
		return _spikyKernelNormCoeff1 * x * x;
	}
	else return 0;
}

template class SmoothedParticles<2>;
template class SmoothedParticles<3>;

}
