#include "SmthParticleHydrodLiquid.h"

#include "Utilities/Constants.h"

namespace PhysX {

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::writeDescription(YAML::Node &root) const
{
	{ // Description of particles.
		YAML::Node node;
		node["name"] = "particles";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["material"]["diffuse_albedo"] = Vector4f(0, 0, 0, 1);
		node["indexed"] = false;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::advance(const real dt)
{
	moveParticles(dt);

	_particles.resetNearbySearcher();
	_particles.computeDensities();

	applyPressureForce(dt);
	applyExternalForces(dt);
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::moveParticles(const real dt)
{
	_particles.parallelForEach([&](const int i) {
		_particles.positions[i] += _velocities[i] * dt;
	});

	// Resolve collisions.
	for (const auto &collider : _colliders) {
		collider->collide(_particles, _velocities);
	}
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::applyPressureForce(const real dt)
{
	// Compute pressures by the equation of state.
	const real eosScale = _targetDensity * _speedOfSound * _speedOfSound / _eosExponent;
	_particles.parallelForEach([&](const int i) {
		_pressures[i] = eosScale * (std::pow(_particles.densities[i] / _targetDensity, _eosExponent) - 1);
		if (_pressures[i] < 0) _pressures[i] = 0;
	});

	// Apply pressure gradient.
	_particles.parallelForEach([&](const int i) {
		_velocities[i] -= _pressures.gradientAtDataPoint(i) / _particles.densities[i] * dt;
	});
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::applyExternalForces(const real dt)
{
	if (_enableGravity) {
		_particles.parallelForEach([&](const int i) {
			_velocities[i][1] -= kGravity * dt;
		});
	}
}

template class SmthParticleHydrodLiquid<2>;
template class SmthParticleHydrodLiquid<3>;

}
