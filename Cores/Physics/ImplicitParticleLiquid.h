#pragma once

#include "Physics/ParticleInCellLiquid.h"
#include "Structures/ParticlesAttribute.h"

namespace PhysX {

template <int Dim>
class ImplicitParticleLiquid : public ParticleInCellLiquid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using EulerianFluid<Dim>::_grid;
	using EulerianFluid<Dim>::_velocity;
	using ParticleInCellLiquid<Dim>::_markerPositions;
	using ParticleInCellLiquid<Dim>::_markerVelocities;

	const real _propOfPic;
	StaggeredGridBasedVectorField<Dim> _deltaVelocity;

public:

	ImplicitParticleLiquid(const StaggeredGrid<Dim> &grid, const int markersCntPerSubcell, const real propOfPic);

	ImplicitParticleLiquid(const ImplicitParticleLiquid &rhs) = delete;
	ImplicitParticleLiquid &operator=(const ImplicitParticleLiquid &rhs) = delete;
	virtual ~ImplicitParticleLiquid() = default;

	virtual void saveFrame(const std::string &frameDir) const override;
	virtual void loadFrame(const std::string &frameDir) override;

protected:

	virtual void transferFromGridsToParticles() override;
	virtual void transferFromParticlesToGrids() override;
};

}
