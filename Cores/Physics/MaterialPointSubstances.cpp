#include "MaterialPointSubstances.h"

#include "Utilities/Constants.h"

namespace PhysX {

template <int Dim>
void MaterialPointSubstances<Dim>::advance(const real dt)
{
	transferFromGridToParticles(dt);
	advect(dt);
	applyLagrangianForces(dt);
	transferFromParticlesToGrid(dt);
	applyEulerianForces(dt);
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
void MaterialPointSubstances<Dim>::applyLagrangianForces(const real dt)
{
	for (auto &substance : _substances) {
		substance.particles.parallelForEach([&](const int i) {
			const VectorDr &pos = substance.particles.positions[i];
			const MatrixDr &velDrv = substance.velocityDerivatives[i];
			MatrixDr &defGrad = substance.deformationGradients[i];
			VectorDr &vel = substance.velocities[i];

			defGrad = (MatrixDr::Identity() + velDrv * dt) * defGrad;
			Eigen::JacobiSVD<MatrixDr> svd(defGrad, Eigen::ComputeFullU | Eigen::ComputeFullV);
			VectorDr svdS = svd.singularValues();
			MatrixDr svdU = svd.matrixU();
			MatrixDr svdV = svd.matrixV();

			const real jacobian = svdS.prod();

			// Currently only support liquid.
			defGrad = MatrixDr::Identity() * std::sqrt(jacobian);
			const MatrixDr stress = 2 * mu * (defGrad - svdU * svdV.transpose()) * defGrad.transpose() + MatrixDr::Identity() * lambda * jacobian * (jacobian - 1);
			vel -= stress * dt * vol * 4 * _grid.invSpacing() * _grid.invSpacing();
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::applyEulerianForces(const real dt)
{
	if (_enableGravity) {
		_velocity.parallelForEach([&](const VectorDi &node) {
			_velocity[node][1] -= kGravity * dt;
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromGridToParticles(const real dt)
{
	for (auto &substance : _substances) {
		// Clear particle attributes.
		substance.velocities.setZero();
		substance.velocityDerivatives.setZero();

		substance.particles.parallelForEach([&](const int i) {
			VectorDr &pos = substance.particles.positions[i];
			VectorDr &vel = substance.velocities[i];
			MatrixDr &velDrv = substance.velocityDerivatives[i];
			MatrixDr &defmGrad = substance.deformationGradients[i];
			// Transfer into velocities and velocity derivatives.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				vel += weight * _velocity[node];
				velDrv += _velocity[node] * deltaPos.transpose() * 4 * _velocity.invSpacing() * weight;
			}
			// Advect particles.
			pos += vel * dt;
			defmGrad = (MatrixDr::Identity() + velDrv * dt) * defmGrad;
		});
	}
}

template <int Dim>
void MaterialPointSubstances<Dim>::transferFromParticlesToGrid(const real dt)
{
	_velocity.setZero();

	for (auto &substance : _substances) {
		const real mass = substance.particles.mass();
		const real stressCoeff = dt * 4 * _velocity.invSpacing() * _velocity.invSpacing() * mass / substance.density;

		substance.particles.forEach([&](const int i) {
			const VectorDr pos = substance.particles.positions[i];
			const VectorDr vel = substance.velocities[i];
			const MatrixDr velDrv = substance.velocityDerivatives[i];
			const MatrixDr stress = substance.computeCauchyStressTensor(i) * stressCoeff;
			// Transfer into velocity.
			for (const auto [node, weight] : _velocity.grid()->quadraticBasisSplineIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity.position(node) - pos;
				_velocity[node] += (vel * mass + (velDrv * mass + stress) * deltaPos) * weight;
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
