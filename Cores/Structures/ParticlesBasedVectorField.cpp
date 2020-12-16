#include "ParticlesBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> ParticlesBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr val = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		val += _data[j] * particles->kernel(nearbyPos - pos) / particles->densities[j];
	});
	return val * particles->mass();
}

template <int Dim>
real ParticlesBasedVectorField<Dim>::divergence(const VectorDr &pos) const
{
	real div = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		div += _data[j].dot(particles->gradientKernel(nearbyPos - pos)) / particles->densities[j];
	});
	return div * particles->mass();
}

template <int Dim>
real ParticlesBasedVectorField<Dim>::differenceDivergenceAtDataPoint(const int idx) const
{
	real div = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		div += (_data[j] - _data[idx]).dot(particles->gradientKernel(nearbyPos - pos)) / particles->densities[j];
	});
	return div * particles->mass();
}

template <int Dim>
real ParticlesBasedVectorField<Dim>::symmetricDivergenceAtDataPoint(const int idx) const
{
	real div = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	const VectorDr posValDivBySquaredDensity = _data[idx] / (particles->densities[idx] * particles->densities[idx]);
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		div += (posValDivBySquaredDensity + _data[j] / (particles->densities[j] * particles->densities[j])).dot(particles->gradientKernel(nearbyPos - pos));
	});
	return div * particles->mass() * particles->densities[idx];
}

template<int Dim>
Vector<Dim, real> ParticlesBasedVectorField<Dim>::laplacianAtDataPoint(const int idx) const
{
	VectorDr lapl = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		lapl += (_data[j] - _data[idx]) * particles->improvedLaplacianKernel(nearbyPos - pos) / particles->densities[j];
	});
	return lapl * particles->mass();
}

template<int Dim>
Vector<Dim, real> ParticlesBasedVectorField<Dim>::divFreeLaplacianAtDataPoint(const int idx) const
{
	VectorDr lapl = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	const VectorDr pos = particles->positions[idx];
	particles->forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
		const VectorDr deltaPos = nearbyPos - pos;
		if (deltaPos.any()) {
			lapl += (_data[j] - _data[idx]).dot(deltaPos) / deltaPos.squaredNorm() * particles->gradientKernel(deltaPos) / particles->densities[j];
		}
	});
	return lapl * 2 * (Dim + real(2)) * particles->mass();
}

template class ParticlesBasedVectorField<2>;
template class ParticlesBasedVectorField<3>;

}
