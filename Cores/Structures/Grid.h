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

	VectorDr cellOrigin() const { return _origin + VectorDr::Ones() * real(0.5) * _spacing; }
	VectorDr faceOrigin(const int axis) const { return _origin + (VectorDr::Ones() - VectorDr::Unit(axis)) * real(0.5) * _spacing; }

	VectorDi cellSize() const { return _resolution; }
	VectorDi faceSize(const int axis) const { return _resolution + VectorDi::Unit(axis); }

	VectorDr cellCenter(const VectorDi &coord) const { return cellOrigin() + coord.cast<real>() * _spacing; }
	VectorDr faceCenter(const int axis, const VectorDi &coord) const { return faceOrigin(axis) + coord.cast<real>() * _spacing; }

protected:

	void getCellLerp(const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const { getLerp(cellOrigin(), pos, coords, weights); }
	void getFaceLerp(const int axis, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const { getLerp(faceOrigin(axis), pos, coords, weights); }

private:

	void getLerp(const VectorDr &dataOrigin, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const;
};

}
