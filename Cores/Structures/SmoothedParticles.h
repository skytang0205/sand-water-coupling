#pragma once

#include "Structures/Particles.h"

namespace PhysX {

template <int Dim>
class SmoothedParticles final : public Particles<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _radius;

public:

	SmoothedParticles(const real mass, const real radius, const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero()) :
		Particles<Dim>(mass, cnt, pos),
		_radius(radius)
	{ }

	SmoothedParticles &operator=(const SmoothedParticles &rhs) = delete;
	virtual ~SmoothedParticles() = default;

	real radius() const { return _radius; }
};

}
