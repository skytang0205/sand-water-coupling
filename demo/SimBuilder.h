#pragma once

#include "Simulation.h"

namespace Pivot {
	struct SimBuildOptions {
		Simulation::Scene Scene;
		int               Scale;
		double            ParticleRadius;
	};

	class SimBuilder {
	public:
		static std::unique_ptr<Simulation> Build(SimBuildOptions const &options);
	
	private:
		static std::unique_ptr<Simulation> BuildFalling(SimBuildOptions const &options);
		static std::unique_ptr<Simulation> BuildBigBall(SimBuildOptions const &options);
		static std::unique_ptr<Simulation> BuildSlope  (SimBuildOptions const &options);

		static void AddParticles(Simulation *simulation, Surface const &surface, bool if_Poission = false, std::function<Vector3d(Vector3d const &)> velocity = nullptr);
	};
}
