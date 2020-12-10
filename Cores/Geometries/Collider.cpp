#include "Collider.h"

namespace PhysX {

template <int Dim>
void Collider<Dim>::collide(Particles<Dim> &particles, ParticlesBasedVectorData<Dim> &velocities, const real radius)
{
	particles.parallelForEach([&](const int i) {
		collide(particles.positions[i], velocities[i], radius);
	});
}

template<int Dim>
void Collider<Dim>::collide(VectorDr &pos, VectorDr &vel, const real radius)
{
	if (surface()->isInside(surface()->signedDistance(pos) - radius)) {
		const VectorDr normal = surface()->closestNormal(pos);
		pos = surface()->closestPosition(pos) + radius * normal; // geometric fix

		const VectorDr colliderVel = velocityAt(pos);
		const VectorDr relativeVel = vel - colliderVel;
		const real normalDotRelativeVel = normal.dot(relativeVel);
		// Check if the velocity is facing opposite direction of the surface normal.
		if (normalDotRelativeVel < 0) {
			// Apply restitution coefficient to the surface normal component of the velocity.
			VectorDr relativeVelN = normalDotRelativeVel * normal;
			VectorDr relativeVelT = relativeVel - relativeVelN;
			const VectorDr deltaRelativeVelN = -(_restitutionCoefficient + 1) * relativeVelN;
			relativeVelN *= -_restitutionCoefficient;
			// Apply friction to the tangential component of the velocity.
			if (relativeVelT.any()) {
				relativeVelT *= std::max(1 - _frictionCoefficient * deltaRelativeVelN.norm() / relativeVelT.norm(), real(0));
			}
			// Reassemble the components.
			vel = relativeVelN + relativeVelT + colliderVel;
		}
	}
}

template class Collider<2>;
template class Collider<3>;

template class StaticCollider<2>;
template class StaticCollider<3>;

template class DynamicCollider<2>;
template class DynamicCollider<3>;

}
