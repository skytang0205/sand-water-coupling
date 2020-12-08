#pragma once

#include "Physics/Simulation.h"

namespace PhysX {

template <int Dim>
class SmthParticleHydrodLiquid : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

public:

	SmthParticleHydrodLiquid(const SmthParticleHydrodLiquid &rhs) = delete;
	SmthParticleHydrodLiquid &operator=(const SmthParticleHydrodLiquid &rhs) = delete;
	virtual ~SmthParticleHydrodLiquid() = default;
};

}
