#pragma once

#include "Utilities/Types.h"

#include <array>
#include <functional>

namespace PhysX {

inline constexpr int CellCentered = 0;
inline constexpr int FaceCentered = 1;
inline constexpr int VertexCentered = 2;

template <int Dim>
class Grid
{
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

	virtual ~Grid() = default;

	real spacing() const { return _spacing; }

	VectorDr cellOrigin() const { return _origin + VectorDr::Ones() * real(0.5) * _spacing; }
	VectorDr faceOrigin(const int axis) const { return _origin + (VectorDr::Ones() - VectorDr::Unit(axis)) * real(0.5) * _spacing; }

	VectorDi cellSize() const { return _resolution; }
	VectorDi faceSize(const int axis) const { return _resolution + VectorDi::Unit(axis); }

	size_t cellCnt() const { return cellSize().cast<size_t>().prod(); }
	size_t faceCnt(const int axis) const { return faceSize(axis).cast<size_t>().prod(); }

	VectorDr cellCenter(const VectorDi &cell) const { return cellOrigin() + cell.cast<real>() * _spacing; }
	VectorDr faceCenter(const int axis, const VectorDi &face) const { return faceOrigin(axis) + face.cast<real>() * _spacing; }

	void getCellLerp(const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const { getLerp(cellOrigin(), cellSize(), pos, coords, weights); }
	void getFaceLerp(const int axis, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const { getLerp(faceOrigin(axis), faceSize(axis), pos, coords, weights); }

	void forEachCell(const std::function<void(const VectorDi &)> &func) const { forEach(cellSize(), func); }
	void parallelForEachCell(const std::function<void(const VectorDi &)> &func) const { parallelForEach(cellSize(), func); };
	void forEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) forEach(faceSize(axis), std::bind(func, axis, std::placeholders::_1)); }
	void parallelForEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) parallelForEach(faceSize(axis), std::bind(func, axis, std::placeholders::_1)); }

protected:

	void getLerp(const VectorDr &dataOrigin, const VectorDi &dataSize, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const;

	void forEach(const VectorDi &dataSize, const std::function<void(const VectorDi &)> &func) const;
	void parallelForEach(const VectorDi &dataSize, const std::function<void(const VectorDi &)> &func) const;
};

}
