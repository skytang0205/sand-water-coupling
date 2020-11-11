#pragma once

#include "Utilities/Types.h"

#include <array>

namespace PhysX {

inline constexpr int CellCentered = 0;
inline constexpr int FaceCentered = 1;
inline constexpr int VertexCentered = 2;

template <int Dim>
class Grid
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _spacing;
	const VectorDi _resolution;
	const VectorDr _origin;

public:

	Grid(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero()) :
		_spacing(spacing),
		_resolution(resolution),
		_origin(center - resolution.cast<real>() * spacing * real(0.5))
	{ }

	Grid(const Grid &rhs) = default;
	Grid &operator=(const Grid &rhs) = default;
	virtual ~Grid() = default;

	VectorDr cellCenter(const VectorDi &coord) const { return _origin + (coord.cast<real>() + VectorDr::Ones() * real(0.5)) * _spacing; }
	VectorDr faceCenter(const int axis, const VectorDi &coord) const { return _origin + (coord.cast<real>() + (VectorDr::Ones() - VectorDr::Unit(axis)) * real(0.5)) * _spacing; }

protected:

	VectorDi getLowerCell(const VectorDr &pos, VectorDr &frac) const
	{
		const VectorDi cell = ((pos - _origin) / _spacing - VectorDr::Ones() * real(0.5)).cast<int>();
		frac = (pos - _origin - (cell.cast<real>() + VectorDr::Ones() * real(0.5)) * _spacing) / _spacing;
		return cell;
	}
};

}
