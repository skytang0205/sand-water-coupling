#pragma once

#include "Utilities/Types.h"
#include "Utilities/IO.h"

#include <functional>
#include <vector>

namespace PhysX {

template <int Dim>
class Particles
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::vector<VectorDr> _positions;

public:

	Particles(const size_t cnt = 0, const VectorDr &pos = VectorDr::Zero()) { reset(cnt, pos); }

	virtual ~Particles() = default;

	virtual void reset(const size_t cnt, const VectorDr &pos = VectorDr::Zero()) { _positions.resize(cnt, pos); }

	size_t size() const { return _positions.size(); }
	virtual void clear() { _positions.clear(); }
	bool empty() const { return _data.empty(); }
	virtual void add(const VectorDr &pos = VectorDr::Zero()) { _positions.push_back(pos); }

	VectorDr *positionsData() { return _positions.data(); }
	const VectorDr *positionsData() const { return _positions.data(); }

	VectorDr &position(const size_t index) { return _positions[index]; }
	const VectorDr &position(const size_t index) const { return _positions[index]; }

	void forEach(const std::function<void(const int)> &func) const
	{
		for (int i = 0; i < _positions.size(); i++) func(i);
	}

	void parallelForEach(const std::function<void(const int)> &func) const
	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < _positions.size(); i++) func(i);
	}

	virtual void load(std::istream &in) { IO::readArray(in, _positions.data(), _positions.size()); }
	virtual void save(std::ostream &out) const { IO::writeArray(out, _positions.data(), _positions.size()); }
};

}
