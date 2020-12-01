#pragma once

#include "Utilities/Types.h"
#include "Utilities/IO.h"

#include <vector>

namespace PhysX {

template <int Dim, typename Type = real>
class ParticlesAttribute
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::vector<Type> _data;

public:

	ParticlesAttribute(const size_t cnt = 0, const Type &val = Zero<Type>()) { resize(cnt, val); }

	virtual ~ParticlesAttribute() = default;

	size_t size() const { return _data.size(); }
	void resize(const size_t cnt, const Type &val = Zero<Type>()) { _data.resize(cnt, val); }
	void clear() { _data.clear(); }
	bool empty() const { return _data.empty(); }
	void add(const Type &val) { _data.push_back(val); }

	auto begin() { return _data.begin(); }
	auto begin() const { return _data.begin(); }
	auto end() { return _data.end(); }
	auto end() const { return _data.end(); }

	Type *data() { return _data.data(); }
	const Type *data() const { return _data.data(); }

	Type &operator[](const int idx) { return _data[idx]; }
	const Type &operator[](const int idx) const { return _data[idx]; }

	void setConstant(const Type &value) { std::fill(_data.begin(), _data.end(), value); }
	void setZero() { setConstant(Zero<Type>()); }

	void forEach(const std::function<void(const int)> &func) const
	{
		for (int i = 0; i < _data.size(); i++) func(i);
	}

	void parallelForEach(const std::function<void(const int)> &func) const
	{
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int i = 0; i < _data.size(); i++) func(i);
	}

	void load(std::istream &in) { IO::readArray(in, _data.data(), _data.size()); }
	void save(std::ostream &out) const { IO::writeArray(out, _data.data(), _data.size()); }
};

template <int Dim> using ParticlesScalarAttribute = ParticlesAttribute<Dim, real>;
template <int Dim> using ParticlesVectorAttribute = ParticlesAttribute<Dim, Vector<Dim, real>>;

}
