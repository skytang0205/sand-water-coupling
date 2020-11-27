#pragma once

#include "Utilities/Types.h"

#include <vector>

namespace PhysX {

template <int Dim, typename Type = real>
class ParticlesAttribute
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::vector<Type> _data;

public:

	ParticlesAttribute() = default;
	virtual ~ParticlesAttribute() = default;

	size_t size() const { return _data.size(); }
	void resize(const size_t cnt) { _data.resize(cnt); }
	void clear() { _data.clear(); }
	void add(const Type &val) { _data.push_back(val); }

	Type &operator[](const size_t idx) { return _data[idx]; }
	const Type &operator[](const size_t idx) const { return _data[idx]; }

	void forEach(const std::function<void(const size_t)> &func) const
	{
		for (size_t i = 0; i < _data.size(); i++) func(i);
	}

	void parallelForEach(const std::function<void(const size_t)> &func) const
	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (size_t i = 0; i < _data.size(); i++) func(i);
	}
};

template <int Dim> using ParticlesScalarAttribute = ParticlesAttribute<Dim, real>;
template <int Dim> using ParticlesVectorAttribute = ParticlesAttribute<Dim, VectorDr>;

}
