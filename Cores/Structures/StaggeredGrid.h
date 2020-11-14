#pragma once

#include "Structures/Grid.h"

namespace PhysX {

template <int Dim>
class StaggeredGrid final
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _spacing;
	const VectorDi _resolution;
	const VectorDr _origin;

	const Grid<Dim> _nodeGrid;
	const Grid<Dim> _cellGrid;
	const std::array<Grid<Dim>, Dim> _faceGrids;

public:

	StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero());

	virtual ~StaggeredGrid() = default;

	real spacing() const { return _spacing; }
	VectorDi resolution() const { return _resolution; }
	VectorDr origin() const { return _origin; }

	const Grid<Dim> *nodeGrid() const { return &_nodeGrid; }
	const Grid<Dim> *cellGrid() const { return &_cellGrid; }
	const Grid<Dim> *faceGrid(const int axis) const { return &_faceGrids[axis]; }

	VectorDi nodeSize() const { return _nodeGrid.dataSize(); }
	VectorDi cellSize() const { return _cellGrid.dataSize(); }
	VectorDi faceSize(const int axis) const { return _faceGrids[axis].dataSize(); }

	size_t nodeCount() const { return _nodeGrid.dataCount(); }
	size_t cellCount() const { return _cellGrid.dataCount(); }
	size_t faceCount(const int axis) const { return _faceGrids[axis].dataCount(); }

	VectorDr nodeCenter(const VectorDi &node) const { return _nodeGrid.dataPosition(node); }
	VectorDr cellCenter(const VectorDi &cell) const { return _cellGrid.dataPosition(cell); }
	VectorDr faceCenter(const int axis, const VectorDi &face) const { return _faceGrids[axis].dataPosition(face); }

	bool isBoundaryFace(const int axis, const VectorDi &face) const { return face[axis] == 0 || face[axis] == _resolution[axis]; }

	void forEachNode(const std::function<void(const VectorDi &)> &func) const { _nodeGrid.forEach(func); }
	void forEachCell(const std::function<void(const VectorDi &)> &func) const { _cellGrid.forEach(func); }
	void forEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) _faceGrids[axis].forEach(std::bind(func, axis, std::placeholders::_1)); }

	void parallelForEachNode(const std::function<void(const VectorDi &)> &func) const { _nodeGrid.parallelForEach(func); }
	void parallelForEachCell(const std::function<void(const VectorDi &)> &func) const { _cellGrid.parallelForEach(func); }
	void parallelForEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) _faceGrids[axis].parallelForEach(std::bind(func, axis, std::placeholders::_1)); }
};

}
