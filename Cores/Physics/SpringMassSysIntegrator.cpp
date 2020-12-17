#include "SpringMassSysIntegrator.h"

namespace PhysX {

template <int Dim>
void SpringMassSysIntegrator<Dim>::accumulateAccelerations(
	const ParticlesVectorAttribute<Dim> &positions,
	const ParticlesVectorAttribute<Dim> &velocities,
	ParticlesVectorAttribute<Dim> &accelerations,
	const ParticlesVectorAttribute<Dim> *const externalForces,
	const std::unordered_set<int> *const constrainedDofs) const
{
	accelerations.setZero();

	for (const auto &spring : *_springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = positions[pid1] - positions[pid0];
		const VectorDr v01 = velocities[pid1] - velocities[pid0];
		const real length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const VectorDr force = e01 * ((length - spring.restLength) * spring.stiffnessCoeff + e01.dot(v01) * spring.dampingCoeff);
		accelerations[pid0] += force * _particles->invMass();
		accelerations[pid1] -= force * _particles->invMass();
	}

	if (externalForces)
		accelerations.asVectorXr() += externalForces->asVectorXr() * _particles->invMass();

	if (constrainedDofs) {
		// Constrain degrees of freedom.
		for (const int dof : *constrainedDofs) {
			*(reinterpret_cast<real *>(accelerations.data()) + dof) = 0;
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
	accumulateAccelerations(positions, velocities, _accelerations, externalForces, constrainedDofs);
	velocities.asVectorXr() += _accelerations.asVectorXr() * dt;
}

template <int Dim>
SmsBackwardEulerIntegrator<Dim>::SmsBackwardEulerIntegrator() :
	_solver(std::make_unique<IcPCgSolver>())
{ }

template <int Dim>
void SmsBackwardEulerIntegrator<Dim>::integrate(
	const ParticlesVectorAttribute<Dim> &positions,
	ParticlesVectorAttribute<Dim> &velocities,
	const real dt,
	const ParticlesVectorAttribute<Dim> *const externalForces,
	const std::unordered_set<int> *const constrainedDofs)
{
	_coeffBackwardEuler.clear();

	const auto addSparseBlockValue = [&](const int pid0, const int pid1, const MatrixDr &block) {
		const int offset0 = pid0 * Dim;
		const int offset1 = pid1 * Dim;
		for (int i = 0; i < Dim; i++) {
			if (constrainedDofs && constrainedDofs->contains(offset0 + i)) continue;
			for (int j = 0; j < Dim; j++) {
				if (constrainedDofs && constrainedDofs->contains(offset1 + j)) continue;
				_coeffBackwardEuler.push_back(Tripletr(offset0 + i, offset1 + j, block(i, j) * _particles->invMass()));
			}
		}
	};

	// Set diagonal blocks.
	_particles->forEach([&](const int i) {
		addSparseBlockValue(i, i, MatrixDr::Identity() * _particles->mass());
	});

	// Set Jacobian of accelaration to velocity.
	for (const auto &spring : *_springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = positions[pid1] - positions[pid0];
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.dampingCoeff * e01 * e01.transpose() * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}

	_matBackwardEuler.setFromTriplets(_coeffBackwardEuler.begin(), _coeffBackwardEuler.end());

	accumulateAccelerations(positions, velocities, _accelerations, externalForces, constrainedDofs);
	_rhsBackwardEuler = _matBackwardEuler * velocities.asVectorXr() + _accelerations.asVectorXr() * dt;

	// Set Jacobian of acceleration to position.
	for (const auto &spring : *_springs) {
		const int pid0 = spring.pid0;
		const int pid1 = spring.pid1;
		const VectorDr r01 = positions[pid1] - positions[pid0];
		const double length = r01.norm();
		const VectorDr e01 = r01.normalized();
		const MatrixDr block = spring.stiffnessCoeff * ((spring.restLength / length - 1) * MatrixDr::Identity() - spring.restLength / length * e01 * e01.transpose()) * dt * dt;
		addSparseBlockValue(pid0, pid0, -block);
		addSparseBlockValue(pid0, pid1, block);
		addSparseBlockValue(pid1, pid0, block);
		addSparseBlockValue(pid1, pid1, -block);
	}

	_matBackwardEuler.setFromTriplets(_coeffBackwardEuler.begin(), _coeffBackwardEuler.end());
	_solver->solve(_matBackwardEuler, velocities.asVectorXr(), _rhsBackwardEuler);
}

template class SpringMassSysIntegrator<2>;
template class SpringMassSysIntegrator<3>;
template class SmsSymplecticEulerIntegrator<2>;
template class SmsSymplecticEulerIntegrator<3>;
template class SmsBackwardEulerIntegrator<2>;
template class SmsBackwardEulerIntegrator<3>;

}
