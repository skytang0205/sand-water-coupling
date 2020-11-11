#pragma once

#include "Utilities/Types.h"

#include <vector>

namespace PhysX {

template <int Dim, typename Type = real>
class Array final
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

protected:

	VectorDi _size;
	std::vector<Type> _data;

public:

	Array(const VectorDi &size = VectorDi::Zero(), const Type &value = Type()) :
		_size(size),
		_data(size.cast<size_t>().prod(), value)
	{ }

	Array(const Array &rhs) = default;
	Array &operator=(const Array &rhs) = default;
	virtual ~Array() = default;

	VectorDi size() const { return _size; }

	void resize(const VectorDi &size, const Type &value = Type())
	{
		_size = size;
		_data.resize(size.cast<size_t>().prod(), value);
	}

	const Type &operator[](const VectorDi &coord) const { return _data[index(coord)]; }
	Type &operator[](const VectorDi &coord) { return _data[index(coord)]; }

	Type lerp(const VectorDi &coord, const VectorDr &frac) const
	{
		if constexpr (Dim == 2) {
			return operator[](coord) * (1 - frac.x()) * (1 - frac.y())
				+ operator[](coord + Vector2i(0, 1)) * (1 - frac.x()) * frac.y()
				+ operator[](coord + Vector2i(1, 0)) * frac.x() * (1 - frac.y())
				+ operator[](coord + Vector2i(1, 1)) * frac.x() * frac.y();
		}
		else {
			return operator[](coord) * (1 - frac.x()) * (1 - frac.y()) * (1 - frac.z())
				+ operator[](coord + Vector3i(0, 0, 1)) * (1 - frac.x()) * (1 - frac.y()) * frac.z()
				+ operator[](coord + Vector3i(0, 1, 0)) * (1 - frac.x()) * frac.y() * (1 - frac.z())
				+ operator[](coord + Vector3i(0, 1, 1)) * (1 - frac.x()) * frac.y() * frac.z()
				+ operator[](coord + Vector3i(1, 0, 0)) * frac.x() * (1 - frac.y()) * (1 - frac.z())
				+ operator[](coord + Vector3i(1, 0, 1)) * frac.x() * (1 - frac.y()) * frac.z()
				+ operator[](coord + Vector3i(1, 1, 0)) * frac.x() * frac.y() * (1 - frac.z())
				+ operator[](coord + Vector3i(1, 1, 1)) * frac.x() * frac.y() * frac.z();
		}
	}

protected:

	size_t index(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return coord.x() + size_t(_size.x()) * coord.y();
		else return coord.x() + _size.x() * (coord.y() + size_t(_size.y()) * coord.z());
	}
};

}
