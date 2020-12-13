#include "PredCorrIncomprSphLiquid.h"

namespace PhysX {

template <int Dim>
void PredCorrIncomprSphLiquid<Dim>::advance(const real dt)
{
	applyPressureForce(dt);
	moveParticles(dt);

	applyExternalForces(dt);
	applyViscosityForce(dt);
}

template <int Dim>
void PredCorrIncomprSphLiquid<Dim>::reinitializeParticlesBasedData()
{
	SmthParticleHydrodLiquid<Dim>::reinitializeParticlesBasedData();
	_predDensities.resize(&_particles);

	_particles.resetNearbySearcher();
	_particles.computeDensities();
}

template <int Dim>
void PredCorrIncomprSphLiquid<Dim>::applyPressureForce(const real dt)
{
	const real delta = computeDelta(dt);

	do {
		// Predict densities.
		_particles.parallelForEach([&](const int i) {
			const VectorDr pos = _particles.positions[i];
			_predDensities[i] = _particles.densities[i];
			_particles.forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
				_predDensities[i] += (_velocities[i] - _velocities[j]).dot(_particles.gradientKernel(nearbyPos - pos)) * _particles.mass() * dt;
			});
			_pressures[i] = delta * (_predDensities[i] - _targetDensity);
			if (_pressures[i] < 0) _predDensities[i] = 0, _pressures[i] = 0;
		});

		// Correct velocities.
		_particles.parallelForEach([&](const int i) {
			_velocities[i] -= _pressures.symmetricGradientAtDataPoint(i) / _particles.densities[i] * dt;
		});
	} while (_predDensities.absoluteMax() / _targetDensity - 1 >= _kPredCorrErrorRatio);
}

template <int Dim>
real PredCorrIncomprSphLiquid<Dim>::computeDelta(const real dt) const
{
	const auto square = [](const real x)->real { return x * x; };
	real delta = square(_targetDensity / dt / _particles.mass()) / 2;
	if constexpr (Dim == 2) {
		delta /= square(_particles.firstDerivativeKernel(_particles.radius() * 2)) * 6
			+ square(_particles.firstDerivativeKernel(_particles.radius() * 2 * std::numbers::sqrt3)) * 6;
	}
	else {
		delta /= square(_particles.firstDerivativeKernel(_particles.radius() * 2)) * 12
			+ square(_particles.firstDerivativeKernel(_particles.radius() * 2 * std::numbers::sqrt2)) * 6
			+ square(_particles.firstDerivativeKernel(_particles.radius() * 2 * std::numbers::sqrt3)) * 24;
	}
	return delta;
}

template class PredCorrIncomprSphLiquid<2>;
template class PredCorrIncomprSphLiquid<3>;

}
