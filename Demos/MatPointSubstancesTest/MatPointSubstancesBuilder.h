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
		if (scale < 0) scale = 200;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(3, length / scale, resolution);
		auto substances = std::make_unique<MaterialPointSubstances<Dim>>(grid);

		auto liquid = MaterialPointSubstances<Dim>::Substance("liquid", Vector4f(0, 0, 1, 1), 1000, 250, 0);
		substances->sampleParticlesInsideSurface(liquid, ImplicitSphere<Dim>(VectorDr::Zero(), length / 3), nppsc);

		substances->_substances.emplace_back(liquid);
		return substances;
	}

	static void reportError(const std::string &msg)
	{
		std::cerr << fmt::format("Error: [MatPointSubstancesBuilder] encountered {}.\n{}", msg) << std::endl;
		std::exit(-1);
	}
};

}
