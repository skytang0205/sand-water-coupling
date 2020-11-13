#pragma once

#include "Utilities/IO.h"
#include "Utilities/Types.h"

#include <vector>

namespace PhysX {

template <int Dim, typename Type = real>
class Array final
{
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDi _size;
	std::vector<Type> _data;

public:

	Array(const VectorDi &size = VectorDi::Zero(), const Type &value = Type()) :
		_size(size),
		_data(size.cast<size_t>().prod(), value)
	{ }

	virtual ~Array() = default;

	VectorDi size() const { return _size; }

	void resize(const VectorDi &size, const Type &value = Type())
	{
		_size = size;
		_data.resize(size.cast<size_t>().prod(), value);
	}

	const Type &operator[](const VectorDi &coord) const { return _data[index(coord)]; }
	Type &operator[](const VectorDi &coord) { return _data[index(coord)]; }

	void read(std::istream &in) { IO::readArray(in, _data.data(), _data.size()); }
	void write(std::ostream &out) const { IO::writeArray(out, _data.data(), _data.size()); }

protected:

	size_t index(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return coord.x() + size_t(_size.x()) * coord.y();
		else return coord.x() + _size.x() * (coord.y() + size_t(_size.y()) * coord.z());
	}
};

}
