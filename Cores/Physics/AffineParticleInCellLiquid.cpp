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

	_markerVelocities.parallelForEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		for (int axis = 0; axis < Dim; axis++) {
			_markerVelocityDerivatives[axis][i] = VectorDr::Zero();
			for (const auto [face, gradWeight] : _velocity[axis].grid()->gradientLinearIntrplDataPoints(pos)) {
				_markerVelocityDerivatives[axis][i] += _velocity[axis][face] * gradWeight;
			}
		}
	});
}

template <int Dim>
void AffineParticleInCellLiquid<Dim>::transferFromParticlesToGrid(StaggeredGridBasedData<Dim> &weightSum)
{
	_markerVelocities.forEach([&](const int i) {
		const VectorDr pos = _markerPositions[i];
		const VectorDr vel = _markerVelocities[i];
		for (int axis = 0; axis < Dim; axis++) {
			const VectorDr gradVelAxis = _markerVelocityDerivatives[axis][i];
			for (const auto [face, weight] : _velocity[axis].grid()->linearIntrplDataPoints(pos)) {
				const VectorDr dPos = _velocity[axis].position(face) - pos;
				_velocity[axis][face] += (vel[axis] + gradVelAxis.dot(dPos)) * weight;
				weightSum[axis][face] += weight;
			}
		}
	});
}

template <int Dim>
void AffineParticleInCellLiquid<Dim>::reinitializeMarkers()
{
	ParticleInCellLiquid<Dim>::reinitializeMarkers();
	for (int axis = 0; axis < Dim; axis++) {
		_markerVelocityDerivatives[axis].resize(_markerPositions.size());
		_markerVelocityDerivatives[axis].setZero();
	}
}

template class AffineParticleInCellLiquid<2>;
template class AffineParticleInCellLiquid<3>;

}
