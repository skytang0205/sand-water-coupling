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

	ParticlesVectorAttribute<Dim> positions;
	ParticlesScalarAttribute<Dim> masses;

public:

	Particles(const size_t cnt, const VectorDr &pos = VectorDr::Zero()) { reset(cnt, pos); }
	Particles(const size_t cnt, const VectorDr &pos, const real mass) { reset(cnt, pos, mass); }

	Particles() = default;
	virtual ~Particles() = default;

	void reset(const size_t cnt, const VectorDr &pos = VectorDr::Zero()) { positions._data.resize(cnt, pos); }

	void reset(const size_t cnt, const VectorDr &pos, const real mass)
	{
		positions._data.resize(cnt, pos);
		masses._data.resize(cnt, mass);
	}

	void clear()
	{
		positions._data.clear();
		masses._data.clear();
	}

	void add(const VectorDr &pos = VectorDr::Zero()) { positions._data.push_back(pos); }

	void add(const VectorDr &pos, const real mass)
	{
		positions._data.push_back(pos);
		masses._data.push_back(mass);
	}

	size_t size() const { return positions.size(); }
	bool empty() const { return positions.empty(); }

	void forEach(const std::function<void(const int)> &func) const { positions.forEach(func); }
	void parallelForEach(const std::function<void(const int)> &func) const { positions.parallelForEach(func); }
};

}
