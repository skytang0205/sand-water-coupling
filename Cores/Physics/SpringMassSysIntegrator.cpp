#include "SpringMassSysIntegrator.h"

namespace PhysX {

template <int Dim>
void SpringMassSysIntegrator<Dim>::accumulateForces(
	const ParticlesVectorAttribute<Dim> &positions,
	const ParticlesVectorAttribute<Dim> &velocities,
	ParticlesVectorAttribute<Dim> &forces,
	const ParticlesVectorAttribute<Dim> *const externalForces,
	const std::unordered_set<int> *const constrainedDofs) const
{
	forces.setZero();

	for (const auto &spring : *_springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = positions[pid1] - positions[pid0];
		const VectorDr v01 = velocities[pid1] - velocities[pid0];
		const real length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const VectorDr force = e01 * ((length - spring.restLength) * spring.stiffnessCoeff + e01.dot(v01) * spring.dampingCoeff);
		forces[pid0] += force;
		forces[pid1] -= force;
	}

	if (externalForces)
		forces.asVectorXr() += externalForces->asVectorXr();

	if (constrainedDofs) {
		// Constrain degrees of freedom.
		for (const int dof : *constrainedDofs) {
			*(reinterpret_cast<real *>(forces.data()) + dof) = 0;
		}
	}
}

template <int Dim>
void SmsSymplecticEulerIntegrator<Dim>::integrate(
	const ParticlesVectorAttribute<Dim> &positions,
	ParticlesVectorAttribute<Dim> &velocities,
	const real dt,
	const ParticlesVectorAttribute<Dim> *const externalForces,
	const std::unordered_set<int> *const constrainedDofs)
{
	accumulateForces(positions, velocities, _forces, externalForces, constrainedDofs);
	velocities.asVectorXr() += _forces.asVectorXr() * _particles->invMass() * dt;
}

template class SpringMassSysIntegrator<2>;
template class SpringMassSysIntegrator<3>;
template class SmsSymplecticEulerIntegrator<2>;
template class SmsSymplecticEulerIntegrator<3>;

}
