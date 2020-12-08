#pragma once

#include "Structures/ParticlesAttribute.h"

#include <memory>

namespace PhysX {

template <int Dim>
class ParticlesNearbySearcher
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _radius;
	const real _squaredRadius;

public:

	ParticlesNearbySearcher(const real radius) : _radius(radius), _squaredRadius(_radius * _radius) { }

	ParticlesNearbySearcher(const ParticlesNearbySearcher &rhs) = delete;
	ParticlesNearbySearcher &operator=(const ParticlesNearbySearcher &rhs) = delete;
	virtual ~ParticlesNearbySearcher() = default;

	virtual void reset(const ParticlesVectorAttribute<Dim> &positions) { }

	virtual void forEach(const ParticlesVectorAttribute<Dim> &positions, const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func)
	{
		for (int i = 0; i < positions.size(); i++)
			if ((positions[i] - pos).squaredNorm() < _squaredRadius) func(i, positions[i]);
	}
};

template <int Dim>
class HashGridSearcher : public ParticlesNearbySearcher<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

public:

	HashGridSearcher(const real radius) : ParticlesNearbySearcher<Dim>(radius) { }

	HashGridSearcher(const HashGridSearcher &rhs) = delete;
	HashGridSearcher &operator=(const HashGridSearcher &rhs) = delete;
	virtual ~HashGridSearcher() = default;
};

}
