#pragma once

#include "Physics/EulerianFluid.h"

namespace PhysX {

class EulerianFluidBuilder final
{
public:

	template <int Dim>
	static std::unique_ptr<EulerianFluid<Dim>> build(const int scale, const int option = 0)
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
	static std::unique_ptr<EulerianFluid<Dim>> buildCase0(int scale)
	{
		DECLARE_DIM_TYPES(Dim)
		if (scale < 0) scale = 150;
		const real length = real(2);
		const VectorDi resolution = scale * VectorDi::Ones();
		StaggeredGrid<Dim> grid(length / scale, resolution);
		auto fluid = std::make_unique<EulerianFluid<Dim>>(grid);

		fluid->_velocity.parallelForEach([&](const int axis, const VectorDi &face) {
			const VectorDr pos = fluid->_velocity[axis].position(face);
			if constexpr (Dim == 2) fluid->_velocity[axis][face] = (axis == 0 ? -pos.y() : pos.x()) * grid.spacing() * 50;
		});
		fluid->_colliders.push_back(
			std::make_unique<StaticCollider<Dim>>(
				std::make_unique<ImplicitSphere<Dim>>(VectorDr::Zero(), real(0.5))));
		fluid->_domainBoundaryHandler = [=](const int axis, const VectorDi &face)->real {
			return grid.spacing() * (axis == 0 ? face.y() : face.x());
		};

		return fluid;
	}

	static void reportError(const std::string &msg);
};

}
