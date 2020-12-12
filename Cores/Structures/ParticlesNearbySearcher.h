#pragma once

#include "Structures/Grid.h"
#include "Structures/ParticlesAttribute.h"

#include <memory>

namespace PhysX {

template <int Dim>
class ParticlesNearbySearcher
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _kernelRadius;
	const real _squaredKernelRadius;

public:

	ParticlesNearbySearcher(const real kernelRadius) : _kernelRadius(kernelRadius), _squaredKernelRadius(_kernelRadius * _kernelRadius) { }

	ParticlesNearbySearcher(const ParticlesNearbySearcher &rhs) = delete;
	ParticlesNearbySearcher &operator=(const ParticlesNearbySearcher &rhs) = delete;
	virtual ~ParticlesNearbySearcher() = default;

	virtual void reset(const ParticlesVectorAttribute<Dim> &positions) { }

	virtual void forEach(const ParticlesVectorAttribute<Dim> &positions, const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func)
	{
		for (int j = 0; j < positions.size(); j++)
			if ((positions[j] - pos).squaredNorm() < _squaredKernelRadius) func(j, positions[j]);
	}
};

template <int Dim>
class HashGridSearcher : public ParticlesNearbySearcher<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using ParticlesNearbySearcher<Dim>::_kernelRadius;
	using ParticlesNearbySearcher<Dim>::_squaredKernelRadius;

	static constexpr int _kPrime0 = 73856093;
	static constexpr int _kPrime1 = 19349663;
	static constexpr int _kPrime2 = 83492791;

	const Grid<Dim> _grid;

	std::vector<std::vector<int>> _buckets;

public:

	HashGridSearcher(const real kernelRadius) : ParticlesNearbySearcher<Dim>(kernelRadius), _grid(_kernelRadius, VectorDi::Zero(), VectorDr::Zero()) { }

	HashGridSearcher(const HashGridSearcher &rhs) = delete;
	HashGridSearcher &operator=(const HashGridSearcher &rhs) = delete;
	virtual ~HashGridSearcher() = default;

	virtual void reset(const ParticlesVectorAttribute<Dim> &positions) override;

	virtual void forEach(const ParticlesVectorAttribute<Dim> &positions, const VectorDr &pos, const std::function<void(const int, const VectorDr &)> &func) override;

	int getHashKey(const VectorDr &pos) const { return getHashKey((_grid.getQuadraticLower(pos) + VectorDi::Ones()).eval()); }

	int getHashKey(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return int((_kPrime0 * ullong(coord[0]) | _kPrime1 * ullong(coord[1])) % _buckets.size());
		else return int((_kPrime0 * ullong(coord[0]) | _kPrime1 * ullong(coord[1]) | _kPrime2 * ullong(coord[2])) % _buckets.size());
	}
};

}
