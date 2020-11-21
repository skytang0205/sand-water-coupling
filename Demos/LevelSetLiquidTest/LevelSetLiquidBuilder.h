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
		case 0:
			return buildCase0<Dim>(scale);
		case 1:
			return buildCase1<Dim>(scale);
		default:
			reportError("invalid option");
			return nullptr;
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

		ImplicitSphere<Dim> sphere(VectorDr::Zero(), length / 4);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ImplicitPlane<Dim>>(VectorDr::Zero() - VectorDr::Unit(1) * length / 2, VectorDr::Unit(1) - VectorDr::Unit(0))));
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ImplicitPlane<Dim>>(VectorDr::Zero() - VectorDr::Unit(1) * length / 2, VectorDr::Unit(1) + VectorDr::Unit(0))));
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase1(const int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = std::make_unique<LevelSetLiquid<Dim>>(grid);

		ImplicitSphere<Dim> sphere(VectorDr::Zero(), length / 4);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ComplementarySurface<Dim>>(
					std::make_unique<ImplicitSphere<Dim>>(VectorDr::Zero(), length / 2 * 0.9))));
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	static void reportError(const std::string &msg);
};

}
