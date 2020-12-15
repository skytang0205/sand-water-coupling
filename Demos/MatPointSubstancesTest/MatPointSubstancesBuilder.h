#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Physics/MaterialPointSubstances.h"
#include "Structures/StaggeredGrid.h"

#include <fmt/core.h>

#include <memory>

namespace PhysX {

class MatPointSubstancesBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<MaterialPointSubstances<Dim>> build(const int scale, const int option, const int nppsc)
	{
		switch (option) {
		case 0:
			return buildCase0<Dim>(scale, nppsc);
		case 1:
			return buildCase1<Dim>(scale, nppsc);
		default:
			reportError("invalid option");
			return nullptr;
		}
	}

protected:

	template <int Dim>
	static std::unique_ptr<MaterialPointSubstances<Dim>> buildCase0(int scale, const int nppsc)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 64;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(3, length / scale, resolution);
		auto substances = std::make_unique<MaterialPointSubstances<Dim>>(grid);

		auto liquid = MaterialPointSubstances<Dim>::Substance("liquid", Vector4f(6, 133, 135, 255) / 255, real(1e3), real(1e5), 0);
		substances->sampleParticlesInsideSurface(liquid, ImplicitSphere<Dim>(VectorDr::Zero(), length * real(.4)), nppsc);

		substances->_substances.emplace_back(liquid);
		return substances;
	}

	template <int Dim>
	static std::unique_ptr<MaterialPointSubstances<Dim>> buildCase1(int scale, const int nppsc)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 5;
		const real length = real(1);
		const VectorDi resolution = scale * (VectorDi::Ones() * 9 + VectorDi::Unit(0) * 7);
		StaggeredGrid<Dim> grid(3, length / scale / 8, resolution);
		auto substances = std::make_unique<MaterialPointSubstances<Dim>>(grid);

		auto jelly = MaterialPointSubstances<Dim>::Substance("jelly", Vector4f(237, 85, 59, 255) / 255, real(1e3), real(7.5e3), real(1.2e5));
		substances->sampleParticlesInsideSurface(jelly, ImplicitSphere<Dim>(VectorDr::Unit(1) * length * (.15) - VectorDr::Unit(0) * length / 2, length * real(.3)), nppsc);

		auto snow = MaterialPointSubstances<Dim>::Substance("snow", Vector4f(238, 238, 240, 255) / 255, real(1e3), real(2.5e4), real(4e5), real(2.5e-2), real(4.5e-3), 10);
		substances->sampleParticlesInsideSurface(snow, ImplicitSphere<Dim>(VectorDr::Unit(1) * length * (.15) + VectorDr::Unit(0) * length / 2, length * real(.3)), nppsc);

		substances->_substances.emplace_back(jelly);
		substances->_substances.emplace_back(snow);
		return substances;
	}


	static void reportError(const std::string &msg)
	{
		std::cerr << fmt::format("Error: [MatPointSubstancesBuilder] encountered {}.\n{}", msg) << std::endl;
		std::exit(-1);
	}
};

}
