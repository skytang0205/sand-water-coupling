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
		node["material"]["diffuse_albedo"] = (Vector4f(52, 108, 156, 255) / 255).eval(); // Haijun Blue
		node["indexed"] = false;
		root["objects"].push_back(node);
	}
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write particles.
		std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
	}
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::saveFrame(const std::string &frameDir) const
{
	{ // Save particles.
		std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.save(fout);
		_velocities.save(fout);
	}
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::loadFrame(const std::string &frameDir)
{
	{ // Load particles.
		std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.load(fin);
		_velocities.load(fin);
	}
	reinitializeParticlesBasedData();
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::initialize()
{
	reinitializeParticlesBasedData();
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::advance(const real dt)
{
	moveParticles(dt);

	applyExternalForces(dt);
	applyViscosityForce(dt);
	applyPressureForce(dt);

	applyPseudoViscosity(dt);
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::reinitializeParticlesBasedData()
{
	_pressures.resize(&_particles);
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

	_particles.resetNearbySearcher();
	_particles.computeDensities();
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

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::applyViscosityForce(const real dt)
{
	if (_viscosityCoeff) {
		auto newVelocities = _velocities;
		_particles.parallelForEach([&](const int i) {
			newVelocities[i] += _viscosityCoeff * _velocities.laplacianAtDataPoint(i) * dt;
		});
		_velocities = newVelocities;
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

/*	_particles.forEach([&](const int i) {
		std::cout << '(' << _particles.positions[i].transpose() << ") "
			<< _particles.densities[i] << ' ' << _pressures[i] << " ["
			<< _pressures.gradientAtDataPoint(i).transpose() / _particles.densities[i] << ']' << std::endl;
	});*/

	// Apply pressure gradient.
	_particles.parallelForEach([&](const int i) {
		_velocities[i] -= _pressures.gradientAtDataPoint(i) / _particles.densities[i] * dt;
	});
}

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::applyPseudoViscosity(const real dt)
{
	auto smoothedVelocities = _velocities;
	_particles.parallelForEach([&](const int i) {
		const VectorDr pos = _particles.positions[i];
		real weightSum = 0;
		VectorDr smoothedVelocity = VectorDr::Zero();
		_particles.forEachNearby(pos, [&](const int j, const VectorDr &nearbyPos) {
			const real weight = _particles.stdKernel(nearbyPos - pos) / _particles.densities[j];
			weightSum += weight;
			smoothedVelocity += weight * _velocities[j];
		});
		if (weightSum > 0) smoothedVelocity /= weightSum;
		smoothedVelocities[i] = smoothedVelocity;
	});

	const real factor = std::clamp(dt * _pseudoViscosityCoeff, real(0), real(1));
	_particles.parallelForEach([&](const int i) {
		_velocities[i] = (1 - factor) * _velocities[i] + factor * smoothedVelocities[i];
	});
}

template class SmthParticleHydrodLiquid<2>;
template class SmthParticleHydrodLiquid<3>;

}
