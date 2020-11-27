#pragma once

#include "Physics/EulerianFluid.h"
#include "Structures/ParticlesAttribute.h"

namespace PhysX {

template <int Dim>
class ParticleInCellLiquid : public EulerianFluid<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	friend class ParticleInCellLiquidBuilder;

protected:

	ParticlesVectorAttribute<Dim> _markerPositions;
	ParticlesVectorAttribute<Dim> _markerVelocities;

public:

	ParticleInCellLiquid(const StaggeredGrid<Dim> &grid);

	ParticleInCellLiquid(const ParticleInCellLiquid &rhs) = delete;
	ParticleInCellLiquid &operator=(const ParticleInCellLiquid &rhs) = delete;
	virtual ~ParticleInCellLiquid() = default;

	virtual void advance(const real dt) override;

protected:

	virtual void advectFields(const real dt) override;
};

}
