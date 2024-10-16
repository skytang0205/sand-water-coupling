#include "Collider.h"

namespace PhysX {

template <int Dim>
void Collider<Dim>::collide(ParticlesVectorAttribute<Dim> &positions, const real radius) const
{
	positions.parallelForEach([&](const int i) {
		collide(positions[i], radius);
	});
}

template <int Dim>
void Collider<Dim>::collide(ParticlesVectorAttribute<Dim> &positions, ParticlesVectorAttribute<Dim> &velocities, const real radius) const
{
	positions.parallelForEach([&](const int i) {
		collide(positions[i], velocities[i], radius);
	});
}

template <int Dim>
void Collider<Dim>::collide(VectorDr &pos, const real radius) const
{
	if (detect(pos, radius)) {
		const VectorDr normal = surface()->closestNormal(pos);
		pos = surface()->closestPosition(pos) + radius * normal; // geometric fix
	}
}

template<int Dim>
void Collider<Dim>::collide(VectorDr &pos, VectorDr &vel, const real radius) const
{
	if (detect(pos, radius)) {
		const VectorDr normal = surface()->closestNormal(pos);
		resolve(pos, normal, vel);
		pos = surface()->closestPosition(pos) + radius * normal; // geometric fix
	}
}

template <int Dim>
bool Collider<Dim>::resolve(const VectorDr &pos, VectorDr &vel) const
{
	const VectorDr normal = surface()->closestNormal(pos);
	return resolve(pos, normal, vel);
}

template <int Dim>
bool Collider<Dim>::resolve(const VectorDr &pos, const VectorDr &normal, VectorDr &vel) const
{
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
		return true;
	}
	return false;
}

template class Collider<2>;
template class Collider<3>;

template class StaticCollider<2>;
template class StaticCollider<3>;

template class DynamicCollider<2>;
template class DynamicCollider<3>;

}
