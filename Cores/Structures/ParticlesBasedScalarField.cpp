#include "ParticlesBasedScalarField.h"

namespace PhysX {

template <int Dim>
real ParticlesBasedScalarField<Dim>::operator()(const VectorDr &pos) const
{
	real val = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		val += _data[i] * particles->stdKernel(nearbyPos - pos) / particles->densities[i];
	});
	return val * particles->mass();
}

template <int Dim>
Vector<Dim, real> ParticlesBasedScalarField<Dim>::gradient(const VectorDr &pos) const
{
	VectorDr grad = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		grad += _data[i] * particles->gradientSpikyKernel(nearbyPos - pos) / particles->densities[i];
	});
	return grad * particles->mass();
}

template <int Dim>
real ParticlesBasedScalarField<Dim>::laplacian(const VectorDr &pos) const
{
	real lapl = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		lapl += _data[i] * particles->laplacianSpikyKernel(nearbyPos - pos) / particles->densities[i];
	});
	return lapl * particles->mass();
}

template <int Dim>
Vector<Dim, real> ParticlesBasedScalarField<Dim>::gradientAtDataPoint(const int idx) const
{
	VectorDr grad = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	const real posValDivBySquaredDensity = _data[idx] / (particles->densities[idx] * particles->densities[idx]);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		grad += (posValDivBySquaredDensity + _data[i] / (particles->densities[i] * particles->densities[i])) * particles->gradientSpikyKernel(nearbyPos - pos);
	});
	return grad * particles->mass() * particles->densities[idx];
}

template <int Dim>
real ParticlesBasedScalarField<Dim>::laplacianAtDataPoint(const int idx) const
{
	real lapl = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		lapl += (_data[i] - _data[idx]) * particles->laplacianSpikyKernel(nearbyPos - pos) / particles->densities[i];
	});
	return lapl * particles->mass();
}

}
