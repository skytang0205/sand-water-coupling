#pragma once

#include "Physics/Simulation.h"

namespace PhysX {

template <int Dim>
class Projectile : public Simulation
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDr _position;
	VectorDr _velocity;
};

}
