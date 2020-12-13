#include "ParticlesBasedScalarField.h"

namespace PhysX {

template <int Dim>
real ParticlesBasedScalarField<Dim>::operator()(const VectorDr &pos) const
{
	real val = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		val += _data[j] * particles->kernel(nearbyPos - pos) / particles->densities[j];
	});
	return val * particles->mass();
}

template <int Dim>
Vector<Dim, real> ParticlesBasedScalarField<Dim>::gradient(const VectorDr &pos) const
{
	VectorDr grad = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		grad += _data[j] * particles->gradientKernel(nearbyPos - pos) / particles->densities[j];
	});
	return grad * particles->mass();
}

template <int Dim>
real ParticlesBasedScalarField<Dim>::laplacian(const VectorDr &pos) const
{
	real lapl = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		lapl += _data[j] * particles->laplacianKernel(nearbyPos - pos) / particles->densities[j];
	});
	return lapl * particles->mass();
}

template <int Dim>
Vector<Dim, real> ParticlesBasedScalarField<Dim>::differenceGradientAtDataPoint(const int idx) const
{
	VectorDr grad = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		grad += (_data[j] - _data[idx]) * particles->gradientKernel(nearbyPos - pos) / particles->densities[j];
	});
	return grad * particles->mass();
}

template <int Dim>
Vector<Dim, real> ParticlesBasedScalarField<Dim>::symmetricGradientAtDataPoint(const int idx) const
{
	VectorDr grad = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	const real posValDivBySquaredDensity = _data[idx] / (particles->densities[idx] * particles->densities[idx]);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		grad += (posValDivBySquaredDensity + _data[j] / (particles->densities[j] * particles->densities[j])) * particles->gradientKernel(nearbyPos - pos);
	});
	return grad * particles->mass() * particles->densities[idx];
}

template <int Dim>
real ParticlesBasedScalarField<Dim>::laplacianAtDataPoint(const int idx) const
{
	real lapl = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		lapl += (_data[j] - _data[idx]) * particles->improvedLaplacianKernel(nearbyPos - pos) / particles->densities[j];
	});
	return lapl * particles->mass();
}

template class ParticlesBasedScalarField<2>;
template class ParticlesBasedScalarField<3>;

}
