#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Physics/SmthParticleHydrodLiquid.h"
#include "Structures/StaggeredGrid.h"

#include <fmt/core.h>

#include <memory>

namespace PhysX {

class SmthPartHydrodLiquidBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<SmthParticleHydrodLiquid<Dim>> build(const int scale, const int option)
	{
		switch (option) {
		case 0:
			return buildCase0<Dim>(scale);
		default:
			reportError("invalid option");
			return nullptr;
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<SmthParticleHydrodLiquid<Dim>> buildCase0(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 15;
		const real length = real(2);

		const real density = 1000;
		const real liquidLen = length / 2;
		const real spacing = liquidLen / scale;
		const real mass = density * std::pow(liquidLen / scale, Dim);
		auto liquid = std::make_unique<SmthParticleHydrodLiquid<Dim>>(mass, spacing * 2);
		liquid->_targetDensity = density;
		liquid->_eosExponent = 1;
		liquid->_speedOfSound = 5;

		Grid<Dim> grid(spacing, scale * VectorDi::Ones(), VectorDr::Zero() - (liquidLen - spacing) / 2 * VectorDr::Ones());
		grid.forEach([&](const VectorDi &cell) {
			liquid->_particles.add(grid.dataPosition(cell));
		});

		liquid->_velocities.resize(&liquid->_particles);

		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ComplementarySurface<Dim>>(
					std::make_unique<ImplicitBox<Dim>>(-length / 2 * VectorDr::Ones(), length * VectorDr::Ones())),
				0));
		return liquid;
	}

	static void reportError(const std::string &msg)
	{
		std::cerr << fmt::format("Error: [SmthPartHydrodLiquidBuilder] encountered {}.\n{}", msg) << std::endl;
		std::exit(-1);
	}
};

}
