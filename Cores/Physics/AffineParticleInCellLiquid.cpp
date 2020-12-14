#include "AffineParticleInCellLiquid.h"

namespace PhysX {

template <int Dim>
AffineParticleInCellLiquid<Dim>::AffineParticleInCellLiquid(const StaggeredGrid<Dim> &grid, const int markersCntPerSubcell) :
	ParticleInCellLiquid<Dim>(grid, markersCntPerSubcell)
{ }

template <int Dim>
void AffineParticleInCellLiquid<Dim>::transferFromGridToParticles()
{
	ParticleInCellLiquid<Dim>::transferFromGridToParticles();

	_particles.parallelForEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		for (int axis = 0; axis < Dim; axis++) {
			_particleVelocityDerivatives[axis][i] = VectorDr::Zero();
			for (const auto [face, gradWeight] : _velocity[axis].grid()->gradientLinearIntrplDataPoints(pos)) {
				_particleVelocityDerivatives[axis][i] += _velocity[axis][face] * gradWeight;
			}
		}
	});
}

template <int Dim>
void AffineParticleInCellLiquid<Dim>::transferFromParticlesToGrid(StaggeredGridBasedScalarData<Dim> &weightSum)
{
	_particles.forEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		const VectorDr vel = _particleVelocities[i];
		for (int axis = 0; axis < Dim; axis++) {
			const VectorDr gradVelAxis = _particleVelocityDerivatives[axis][i];
			for (const auto [face, weight] : _velocity[axis].grid()->linearIntrplDataPoints(pos)) {
				const VectorDr deltaPos = _velocity[axis].position(face) - pos;
				_velocity[axis][face] += (vel[axis] + gradVelAxis.dot(deltaPos)) * weight;
				weightSum[axis][face] += weight;
			}
		}
	});
}

template <int Dim>
void AffineParticleInCellLiquid<Dim>::reinitializeParticlesBasedData()
{
	ParticleInCellLiquid<Dim>::reinitializeParticlesBasedData();
	for (int axis = 0; axis < Dim; axis++) {
		_particleVelocityDerivatives[axis].resize(&_particles);
		_particleVelocityDerivatives[axis].setZero();
	}
}

template class AffineParticleInCellLiquid<2>;
template class AffineParticleInCellLiquid<3>;

}
