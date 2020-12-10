#pragma once

#include "Physics/AffineParticleInCellLiquid.h"
#include "Physics/FlImplicitParticleLiquid.h"
#include "Physics/ParticleInCellLiquid.h"

namespace PhysX {

class ParticleInCellLiquidBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> build(const int scale, const int option, const int nppsc, const real alpha)
	{
		switch (option) {
		case 0:
			return buildCase0<Dim>(scale, nppsc, alpha);
		case 1:
			return buildCase1<Dim>(scale, nppsc, alpha);
		case 2:
			return buildCase2<Dim>(scale, nppsc, alpha);
		case 3:
			return buildCase3<Dim>(scale, nppsc, alpha);
		case 4:
			return buildCase4<Dim>(scale, nppsc, alpha);
		default:
			reportError("invalid option");
			return nullptr;
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> buildCase0(int scale, int nppsc, const real alpha)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = makeLiquid<Dim>(grid, nppsc, alpha);

		ImplicitSphere<Dim> sphere(VectorDr::Zero(), length / 4);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ImplicitPlane<Dim>>(VectorDr::Zero() - VectorDr::Unit(1) * length / 2, VectorDr::Unit(1) - VectorDr::Unit(0))));
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ImplicitPlane<Dim>>(VectorDr::Zero() - VectorDr::Unit(1) * length / 2, VectorDr::Unit(1) + VectorDr::Unit(0))));
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> buildCase1(int scale, int nppsc, const real alpha)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = makeLiquid<Dim>(grid, nppsc, alpha);

		ImplicitSphere<Dim> sphere(VectorDr::Zero(), length / 4);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ComplementarySurface<Dim>>(
					std::make_unique<ImplicitSphere<Dim>>(VectorDr::Zero(), length / 2))));
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> buildCase2(int scale, int nppsc, const real alpha)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = makeLiquid<Dim>(grid, nppsc, alpha);

		ImplicitPlane<Dim> plane(-VectorDr::Unit(1) * length / 8, VectorDr::Unit(1));
		liquid->_levelSet.unionSurface(plane);
		liquid->_levelSet.intersectSurface(ImplicitBox<Dim>(grid.domainOrigin(), grid.domainLengths()));
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> buildCase3(int scale, int nppsc, const real alpha)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = makeLiquid<Dim>(grid, nppsc, alpha);

		ImplicitPlane<Dim> plane(-VectorDr::Unit(1) * length / 8, VectorDr::Unit(1));
		ImplicitSphere<Dim> sphere(VectorDr::Unit(1) * length / 8, length / 8);
		liquid->_levelSet.unionSurface(plane);
		liquid->_levelSet.unionSurface(sphere);
		liquid->_levelSet.intersectSurface(ImplicitBox<Dim>(grid.domainOrigin(), grid.domainLengths()));
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> buildCase4(int scale, int nppsc, const real alpha)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 90;
		const real length = real(2);
		const VectorDi resolution = scale * (VectorDi::Ones() * 2 + VectorDi::Unit(0));
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto liquid = makeLiquid<Dim>(grid, nppsc, alpha);

		ImplicitPlane<Dim> plane(VectorDr::Unit(0) * length / 2, VectorDr::Unit(0));
		liquid->_levelSet.unionSurface(plane);
		liquid->_levelSet.intersectSurface(ImplicitBox<Dim>(grid.domainOrigin(), grid.domainLengths()));
		return liquid;
	}

	template <int Dim>
	static std::unique_ptr<ParticleInCellLiquid<Dim>> makeLiquid(const StaggeredGrid<Dim> &grid, const int nppsc, real alpha)
	{
		if (alpha >= 1)
			return std::make_unique<ParticleInCellLiquid<Dim>>(grid, nppsc);
		else if (alpha >= 0)
			return std::make_unique<FlImplicitParticleLiquid<Dim>>(grid, nppsc, alpha);
		else
			return std::make_unique<AffineParticleInCellLiquid<Dim>>(grid, nppsc);
	}

	static void reportError(const std::string &msg)
	{
		std::cerr << fmt::format("Error: [ParticleInCellLiquidBuilder] encountered {}.\n{}", msg) << std::endl;
		std::exit(-1);
	}
};

}
