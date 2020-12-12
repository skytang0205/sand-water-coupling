#include "PredCorrIncomprSphLiquid.h"

namespace PhysX {

template <int Dim>
void PredCorrIncomprSphLiquid<Dim>::reinitializeParticlesBasedData()
{
	SmthParticleHydrodLiquid<Dim>::reinitializeParticlesBasedData();
	_predPositions.resize(&_particles);
	_predVelocities.resize(&_particles);
	_densityErrors.resize(&_particles);
}

template <int Dim>
void PredCorrIncomprSphLiquid<Dim>::applyPressureForce(const real dt)
{
	const real delta = computeDelta(dt);
	_pressures.setZero();
	_predVelocities = _velocities;

	for (int iter = 0; iter < _kPredCorrMaxIters; iter++) {
		// Predict positions.
		_particles.parallelForEach([&](const int i) {
			_predPositions[i] = _particles.positions[i] + _predVelocities[i] * dt;
		});

		// Resolve collisions.
		for (const auto &collider : _colliders) {
			collider->collide(_predPositions, _predVelocities, _particles.radius());
		}

		// Compute pressure from density error.
		_particles.parallelForEach([&](const int i) {
			real density = 0;
			_particles.forEachNearby(_particles.positions[i], [&](const int j, const VectorDr &nearbyPos) {
				density += _particles.stdKernel(_predPositions[j] - _predPositions[i]);
			});
			density *= _particles.mass();

			real densityError = density - _targetDensity;
			real pressure = delta * densityError;
			if (pressure < 0) pressure = 0, densityError = 0;

			_particles.densities[i] = density;
			_pressures[i] += pressure;
			_densityErrors[i] = densityError;
		});

		// Predict velocities.
		_particles.parallelForEach([&](const int i) {
			_predVelocities[i] = _velocities[i] - _pressures.gradientAtDataPoint(i) / _particles.densities[i] * dt;
		});

		if (_densityErrors.absoluteMax() < _kPredCorrErrorRatio * _targetDensity) break;
	}

	_velocities = _predVelocities;
}

template <int Dim>
real PredCorrIncomprSphLiquid<Dim>::computeDelta(const real dt) const
{
	const auto square = [](const real x)->real { return x * x; };
	real delta = square(_targetDensity / dt / _particles.mass()) / 2;
	if constexpr (Dim == 2) {
		delta /= square(_particles.firstDerivativeSpikyKernel(_particles.radius() * 2)) * 6
			+ square(_particles.firstDerivativeSpikyKernel(_particles.radius() * 2 * std::numbers::sqrt3)) * 6;
	}
	else {
		delta /= square(_particles.firstDerivativeSpikyKernel(_particles.radius() * 2)) * 12
			+ square(_particles.firstDerivativeSpikyKernel(_particles.radius() * 2 * std::numbers::sqrt2)) * 6
			+ square(_particles.firstDerivativeSpikyKernel(_particles.radius() * 2 * std::numbers::sqrt3)) * 24;
	}
	return delta;
}

template class PredCorrIncomprSphLiquid<2>;
template class PredCorrIncomprSphLiquid<3>;

}
