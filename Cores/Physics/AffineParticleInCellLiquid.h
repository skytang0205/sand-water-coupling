#pragma once

#include "Physics/ParticleInCellLiquid.h"

namespace PhysX {

template <int Dim>
class AffineParticleInCellLiquid : public ParticleInCellLiquid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using EulerianFluid<Dim>::_kExtrapMaxSteps;
	using EulerianFluid<Dim>::_velocity;
	using EulerianFluid<Dim>::_boundaryHelper;
	using LevelSetLiquid<Dim>::_levelSet;
	using ParticleInCellLiquid<Dim>::_markerPositions;
	using ParticleInCellLiquid<Dim>::_markerVelocities;

	std::array<ParticlesVectorAttribute<Dim>, Dim> _markerVelocityDerivatives;

public:

	AffineParticleInCellLiquid(const StaggeredGrid<Dim> &grid, const int markersCntPerSubcell);

	AffineParticleInCellLiquid(const AffineParticleInCellLiquid &rhs) = delete;
	AffineParticleInCellLiquid &operator=(const AffineParticleInCellLiquid &rhs) = delete;
	virtual ~AffineParticleInCellLiquid() = default;

protected:

	virtual void transferFromGridToParticles() override;
	virtual void transferFromParticlesToGrid(StaggeredGridBasedData<Dim> &weightSum) override;

	virtual void reinitializeMarkers() override;
};

}
