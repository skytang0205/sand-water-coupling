#include "SimBuilder.h"
#include "VolumeSampler.h"

#include "CSG.h"

namespace Pivot {
	std::unique_ptr<Simulation> SimBuilder::Build(SimBuildOptions const &options) {
		std::unique_ptr<Simulation> simulation;
		switch (options.Scene) {
			case Simulation::Scene::Falling: simulation = BuildFalling(options); break;
			case Simulation::Scene::BigBall: simulation = BuildBigBall(options); break;
			case Simulation::Scene::Slope  : simulation = BuildSlope  (options); break;
		}
		simulation->m_Scene = options.Scene;
		return simulation;
	}

	std::unique_ptr<Simulation> SimBuilder::BuildFalling(SimBuildOptions const &options) {
		constexpr double length = 1.;
		constexpr int bw = 2;
		int const scale = options.Scale < 0 ? 128 : options.Scale;
		double const radius = options.ParticleRadius < 0 ? .5/double(scale) : options.ParticleRadius;
		StaggeredGrid sgrid(bw, length / (scale - bw * 2), Vector3i(1, 1, 1) * scale);
		auto sim = std::make_unique<Simulation>(sgrid,radius);
		CSG::Union(sim->m_LevelSet, ImplicitPlane (-Vector3d::Unit(1) * length * .15, Vector3d::Unit(1)));
		CSG::Union(sim->m_LevelSet, ImplicitSphere( Vector3d::Unit(1) * length * .05, length * .1));
		return sim;
	}

	std::unique_ptr<Simulation> SimBuilder::BuildBigBall(SimBuildOptions const &options) {
		constexpr double length = 1.;
		constexpr int bw = 2;
		int const scale = options.Scale < 0 ? 128 : options.Scale;
		double const radius = options.ParticleRadius < 0 ? .5/double(scale) : options.ParticleRadius;
		StaggeredGrid sgrid(bw, length / (scale - bw * 2), Vector3i(1, 1, 1) * scale);
		auto sim = std::make_unique<Simulation>(sgrid, radius);
		CSG::Union(sim->m_LevelSet, ImplicitSphere(Vector3d::Zero(), length * .25));
		AddParticles(sim.get(), ImplicitSphere(Vector3d::Zero() * length, .25 * length), true);
		return sim;
	}

	std::unique_ptr<Simulation> SimBuilder::BuildSlope(SimBuildOptions const &options) {
		constexpr double length = 1.;
		constexpr int bw = 2;
		int const scale = options.Scale < 0 ? 128 : options.Scale;
		double const radius = options.ParticleRadius < 0 ? .5/double(scale) : options.ParticleRadius;
		StaggeredGrid sgrid(bw, length / (scale - bw * 2), Vector3i(1, 1, 1) * scale);
		auto sim = std::make_unique<Simulation>(sgrid, radius);
		CSG::Union(sim->m_LevelSet, ImplicitSphere(Vector3d::Zero(), length * .25));
		AddParticles(sim.get(), ImplicitSphere(Vector3d::Zero() * length, .25 * length));
		CSG::Union(sim->m_Collider.LevelSet, ImplicitPlane(Vector3d(-2, -1, 0) * length * .25, Vector3d(1, 4, 0).normalized()));
		return sim;
	}

	void SimBuilder::AddParticles(Simulation *simulation, Surface const &surface, bool if_Poission, std::function<Vector3d(Vector3d const &)> velocity) {
		VolumeSampler sampler(simulation->m_ParticleRadius * 2.);
		auto const positions = sampler.Sample(surface, if_Poission);
		simulation->m_DEMParticles.reserve(simulation->m_DEMParticles.size() + positions.size());
		for (auto const &pos : positions) {
			simulation->m_DEMParticles.push_back({
				.Position = pos,
				.Velocity = velocity ? velocity(pos) : Vector3d::Zero(),
			});
		}
	}
}
