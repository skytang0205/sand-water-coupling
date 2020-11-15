#pragma once

#include "Utilities/Types.h"

#include <array>
#include <functional>

namespace PhysX {

template <int Dim>
class Grid final
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _spacing;
	const VectorDi _dataSize;
	const VectorDr _dataOrigin;

public:

	Grid(const real spacing, const VectorDi &dataSize, const VectorDr &dataOrigin) :
		_spacing(spacing),
		_dataSize(dataSize),
		_dataOrigin(dataOrigin)
	{ }

	virtual ~Grid() = default;

	bool isValid(const VectorDi &coord) const { return (coord.array() >= 0).all() && (coord.array() < _dataSize.array()).all(); }

	real spacing() const { return _spacing; }
	VectorDi dataSize() const { return _dataSize; }
	VectorDr dataOrigin() const { return _dataOrigin; }

	size_t dataCount() const { return _dataSize.cast<size_t>().prod(); }
	VectorDr dataPosition(const VectorDi &coord) const { return _dataOrigin + coord.cast<real>() * _spacing; }

	size_t index(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return coord.x() + size_t(_dataSize.x()) * coord.y();
		else return coord.x() + _dataSize.x() * (coord.y() + size_t(_dataSize.y()) * coord.z());
	}

	void getLerpCoordsAndWeights(const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const;

	void forEach(const std::function<void(const VectorDi &)> &func) const;
	void parallelForEach(const std::function<void(const VectorDi &)> &func) const;
};

}
