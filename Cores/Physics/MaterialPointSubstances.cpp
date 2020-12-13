#include "MaterialPointSubstances.h"

namespace PhysX {

template <int Dim>
void MaterialPointSubstances<Dim>::advance(const real dt)
{
	transferFromGridToParticles();
	advect(dt);
	applyParticleForces(dt);
	transferFromParticlesToGrid();
	applyGridForces(dt);
}

template <int Dim>
void MaterialPointSubstances<Dim>::advect(const real dt)
{
	for (auto &substance : _substances) {
		substance.particles.parallelForEach([&](const int i) {
			substance.particles.positions[i] += substance.velocities[i] * dt;
		});
		// Resolve collisions.
		for (const auto &collider : _colliders) {
			collider->collide(substance.particles.positions, substance.velocities);
		}
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromGridToParticles()
{
	for (auto &substance : _substances) {
		// Clear particle attributes.
		substance.velocities.setZero();
		substance.velocityDerivatives.setZero();

		substance.particles.parallelForEach([&](const int i) {
			VectorDr &pos = substance.particles.positions[i];
			VectorDr &vel = substance.velocities[i];
			MatrixDr &velDrv = substance.velocityDerivatives[i];
			// Transfer into velocities and velocity derivatives.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				vel += weight * _velocity[node];
				velDrv += _velocity[node] * deltaPos.transpose() * 4 * _velocity.invSpacing() * weight;
			}
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromParticlesToGrid()
{
	_velocity.setZero();

	for (auto &substance : _substances) {
		substance.particles.forEach([&](const int i) {
			const VectorDr pos = substance.particles.positions[i];
			const VectorDr vel = substance.velocities[i];
			const MatrixDr velDrv = substance.velocityDerivatives[i];
			// Transfer into velocity.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				_velocity[node] += (vel + velDrv * deltaPos) * substance.particles.mass() * weight;
				_mass[node] += substance.particles.mass() * weight;
			}
		});
	}

	_velocity.parallelForEach([&](const VectorDi &node) {
		if (_mass[node])
			_velocity[node] /= _mass[node];
	});
}

template class MaterialPointSubstances<2>;
template class MaterialPointSubstances<3>;
}
