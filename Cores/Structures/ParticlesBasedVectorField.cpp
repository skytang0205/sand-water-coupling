#include "ParticlesBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> ParticlesBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr val = VectorDr::Zero();
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		val += _data[i] * particles->stdKernel(nearbyPos - pos) / particles->densities[i];
	});
	return val * particles->mass();
}

template <int Dim>
real ParticlesBasedVectorField<Dim>::divergence(const VectorDr &pos) const
{
	real div = 0;
	auto particles = static_cast<const SmoothedParticles<Dim> *>(_particles);
	particles->forEachNearby(pos, [&](const int i, const VectorDr &nearbyPos) {
		div += _data[i].dot(particles->gradientSpikyKernel(nearbyPos - pos)) / particles->densities[i];
	});
	return div * particles->mass();
}

template class ParticlesBasedVectorField<2>;
template class ParticlesBasedVectorField<3>;

}
