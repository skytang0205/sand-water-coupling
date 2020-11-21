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
		case 2:
			return buildCase2<Dim>(scale);
		case 3:
			return buildCase3<Dim>(scale);
		default:
			reportError("invalid option");
			return nullptr;
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase0(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
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
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase1(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = std::make_unique<LevelSetLiquid<Dim>>(grid);

		ImplicitSphere<Dim> sphere(VectorDr::Zero(), length / 4);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ComplementarySurface<Dim>>(
					std::make_unique<ImplicitSphere<Dim>>(VectorDr::Zero(), length / 2))));
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase2(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = std::make_unique<LevelSetLiquid<Dim>>(grid);

		ImplicitPlane<Dim> plane(-VectorDr::Unit(1) * length / 8, VectorDr::Unit(1));
		ImplicitSphere<Dim> sphere(VectorDr::Unit(1) * length / 8, length / 8);
		liquid->_levelSet.unionSurface(plane);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<LevelSetLiquid<Dim>> buildCase3(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 90;
		const real length = real(2);
		const VectorDi resolution = scale * (VectorDi::Ones() * 2 + VectorDi::Unit(0));
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = std::make_unique<LevelSetLiquid<Dim>>(grid);

		ImplicitBox<Dim> box(-VectorDr::Unit(0) * length / 2, VectorDr::Ones() * length);
		liquid->_levelSet.unionSurface(box);
		liquid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real { return 0; };
		return liquid;
	}

	static void reportError(const std::string &msg);
};

}
