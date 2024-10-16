#include "SpringMassSystem.h"

#include "Utilities/Constants.h"

namespace PhysX {

template <int Dim>
SpringMassSystem<Dim>::SpringMassSystem() :
	_integrator(std::make_unique<SmsSemiImplicitIntegrator<Dim>>())
{ }

template <int Dim>
void SpringMassSystem<Dim>::writeDescription(YAML::Node &root) const
{
	{ // Description of particles.
		YAML::Node node;
		node["name"] = "particles";
		node["data_mode"] = "dynamic";
		node["primitive_type"] = "point_list";
		node["indexed"] = false;
		node["material"]["diffuse_albedo"] = Vector4f(1, 0, 0, 1);
		root["objects"].push_back(node);
	}
	{ // Description of springs.
		YAML::Node node;
		node["name"] = "springs";
		node["data_mode"] = "semi-dynamic";
		node["primitive_type"] = "line_list";
		node["indexed"] = true;
		node["material"]["diffuse_albedo"] = Vector4f(0, 0, 1, 1);
		root["objects"].push_back(node);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::writeFrame(const std::string &frameDir, const bool staticDraw) const
{
	{ // Write particles.
		std::ofstream fout(frameDir + "/particles.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			_particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
	}
	{ // Write springs.
		std::ofstream fout(frameDir + "/springs.mesh", std::ios::binary);
		IO::writeValue(fout, uint(_particles.size()));
		_particles.forEach([&](const int i) {
			IO::writeValue(fout, _particles.positions[i].template cast<float>().eval());
		});
		if constexpr (Dim == 3) {
			_particles.forEach([&](const int i) {
				IO::writeValue(fout, VectorDf::Unit(2).eval());
			});
		}
		if (staticDraw) {
			IO::writeValue(fout, 2 * uint(_springs.size()));
			for (const auto &spring : _springs) {
				IO::writeValue(fout, uint(spring.pid0));
				IO::writeValue(fout, uint(spring.pid1));
			}
		}
	}
}

template <int Dim>
void SpringMassSystem<Dim>::saveFrame(const std::string &frameDir) const
{
	{ // Save particles.
		std::ofstream fout(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.save(fout);
	}
	{ // Save velocities.
		std::ofstream fout(frameDir + "/velocities.sav", std::ios::binary);
		_velocities.save(fout);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::loadFrame(const std::string &frameDir)
{
	{ // Load particles.
		std::ifstream fin(frameDir + "/particles.sav", std::ios::binary);
		_particles.positions.load(fin);
	}
	reinitializeParticlesBasedData();
	{ // Load velocities.
		std::ifstream fin(frameDir + "/velocities.sav", std::ios::binary);
		_velocities.load(fin);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::initialize()
{
	reinitializeParticlesBasedData();
	_integrator->reset(&_particles, &_springs);
}

template <int Dim>
void SpringMassSystem<Dim>::advance(const real dt)
{
	moveParticles(dt);

	applyExternalForces(dt);
	applyElasticForce(dt);
}

template <int Dim>
void SpringMassSystem<Dim>::reinitializeParticlesBasedData()
{
	_externalForces.resize(&_particles);
}

template <int Dim>
void SpringMassSystem<Dim>::moveParticles(const real dt)
{
	_particles.positions.asVectorXr() += _velocities.asVectorXr() * dt;

	// Resolve collisions.
	for (const auto &collider : _colliders) {
		collider->collide(_particles.positions, _velocities);
	}
}

template <int Dim>
void SpringMassSystem<Dim>::applyExternalForces(const real dt)
{
	_particles.parallelForEach([&](const int i) {
		if (_enableGravity && !_constrainedDofs.contains(i * Dim + 1))
			_velocities[i][1] -= kGravity * dt;
	});
}

template <int Dim>
void SpringMassSystem<Dim>::applyElasticForce(const real dt)
{
	_integrator->integrate(_particles.positions, _velocities, dt, _constrainedDofs);
}

template class SpringMassSystem<2>;
template class SpringMassSystem<3>;

}
