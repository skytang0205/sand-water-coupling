#pragma once

#include "Structures/ParticlesAttribute.h"

#include <functional>
#include <vector>

namespace PhysX {

template <int Dim>
class Particles
{
	DECLARE_DIM_TYPES(Dim)

public:

	const real _mass;
	const real _invMass;
	ParticlesVectorAttribute<Dim> positions;

public:

	Particles(const real mass = 1, const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero()) : _mass(mass), _invMass(1 / _mass) { resize(cnt, pos); }

	Particles &operator=(const Particles &rhs) = delete;
	virtual ~Particles() = default;

	void resize(const size_t cnt, const VectorDr &pos = VectorDr::Zero()) { positions._data.resize(cnt, pos); }

	real mass() const { return _mass; }
	real invMass() const { return _invMass; }
	void add(const VectorDr &pos = VectorDr::Zero()) { positions._data.push_back(pos); }
	void clear() { positions._data.clear(); }
	size_t size() const { return positions.size(); }
	bool empty() const { return positions.empty(); }

	void forEach(const std::function<void(const int)> &func) const { positions.forEach(func); }
	void parallelForEach(const std::function<void(const int)> &func) const { positions.parallelForEach(func); }
};

}
