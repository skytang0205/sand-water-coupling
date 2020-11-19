#pragma once

#include "Physics/LevelSetLiquid.h"

namespace PhysX {

class LevelSetLiquidBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> build(const int scale, const int option = 0)
	{
		switch (option) {
		case 0: return buildCase0<Dim>(scale);
		default: reportError("invalid option");
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase0(const int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = std::make_unique<LevelSetLiquid<Dim>>(grid);

		ImplicitSphere<Dim> sphere(VectorDr::Zero() - VectorDr::Unit(1) * length / 2, length / 2);
		liquid->_levelSet.setFromSurface(sphere);
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	static void reportError(const std::string &msg);
};

}
