#pragma once

#include "Structures/Particles.h"

namespace PhysX {

template <int Dim>
class SmoothedParticles final : public Particles<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::vector<real> _masses;

public:

	SmoothedParticles(const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero(), const real mass = 1) { reset(cnt, pos, mass); }

	virtual void reset(const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero(), const real mass = 1) override
	{
		Particles<Dim>::reset(cnt, pos);
		_masses.resize(cnt, mass);
	}

	virtual void clear() override
	{
		_masses.clear();
		Particles<Dim>::clear();
	}

	void add(const VectorDr &pos, const real mass)
	{
		Particles<Dim>::add(pos);
		_masses.push_back(mass);
	}

	virtual void add(const VectorDr &pos = VectorDr::Zero()) override { add(pos, 1); }

	real *massesData() { return _masses.data(); }
	const real *massesData() const { return _masses.data(); }

	real &mass(const size_t index) { return _masses[index]; }
	const real &mass(const size_t index) const { return _masses[index]; }
};

}
