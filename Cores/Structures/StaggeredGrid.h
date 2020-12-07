#pragma once

#include "Structures/Grid.h"

namespace PhysX {

template <int Dim>
class StaggeredGrid final
{
	DECLARE_DIM_TYPES(Dim)

protected:

	static constexpr int _kBoundaryWidth = 2;

	const real _spacing;
	const real _invSpacing;
	const VectorDi _resolution;
	const VectorDr _origin;

	const Grid<Dim> _nodeGrid;
	const Grid<Dim> _cellGrid;
	const std::array<Grid<Dim>, Dim> _faceGrids;

public:

	StaggeredGrid(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero());

	StaggeredGrid &operator=(const StaggeredGrid &rhs) = delete;
	virtual ~StaggeredGrid() = default;

	real spacing() const { return _spacing; }
	real invSpacing() const { return _invSpacing; }
	VectorDi resolution() const { return _resolution; }
	VectorDr origin() const { return _origin; }

	VectorDr domainOrigin() const { return _origin + VectorDr::Ones() * _kBoundaryWidth * _spacing; }
	VectorDr domainLengths() const { return (_resolution - VectorDi::Ones() * _kBoundaryWidth * 2).template cast<real>() * _spacing; }

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

	bool isInsideFace(const int axis, const VectorDi &face) const { return _faceGrids[axis].isInside(face, _kBoundaryWidth); }
	bool isBoundaryFace(const int axis, const VectorDi &face) const { return face[axis] <= _kBoundaryWidth || face[axis] >= _resolution[axis] - _kBoundaryWidth || !isInsideFace(axis, face); }

	void forEachNode(const std::function<void(const VectorDi &)> &func) const { _nodeGrid.forEach(func); }
	void forEachCell(const std::function<void(const VectorDi &)> &func) const { _cellGrid.forEach(func); }
	void forEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) _faceGrids[axis].forEach(std::bind(func, axis, std::placeholders::_1)); }

	void parallelForEachNode(const std::function<void(const VectorDi &)> &func) const { _nodeGrid.parallelForEach(func); }
	void parallelForEachCell(const std::function<void(const VectorDi &)> &func) const { _cellGrid.parallelForEach(func); }
	void parallelForEachFace(const std::function<void(const int, const VectorDi &)> &func) const { for (int axis = 0; axis < Dim; axis++) _faceGrids[axis].parallelForEach(std::bind(func, axis, std::placeholders::_1)); }

	static constexpr int numberOfCellNodes() { return 1 << Dim; }
	static constexpr int numberOfCellFaces() { return Dim << 1; }
	static constexpr int numberOfCellEdges() { return Dim * (1 << (Dim - 1)); }
	static constexpr int numberOfFaceNodes() { return 1 << (Dim - 1); }

	static VectorDi cellNode(const VectorDi &cell, const int ord)
	{
		if constexpr (Dim == 2) return cell + VectorDi(ord & 1, ord >> 1 & 1);
		else return cell + VectorDi(ord & 1, ord >> 1 & 1, ord >> 2 & 1);
	};

	static VectorDi faceNode(const int axis, const VectorDi &face, const int ord)
	{
		if constexpr (Dim == 2) return face + VectorDi::Unit(axis ^ 1) * (ord & 1);
		else return face + VectorDi::Unit((axis + 1) % 3) * (ord & 1) + VectorDi::Unit((axis + 2) % 3) * (ord >> 1 & 1);
	}

	static constexpr int cellFaceAxis(const int ord) { return ord >> 1; }
	static constexpr int cellFaceSide(const int ord) { return ord & 1 ? 1 : -1; }
	static VectorDi cellFace(const VectorDi &cell, const int ord) { return cell + VectorDi::Unit(cellFaceAxis(ord)) * (ord & 1); }
	static VectorDi faceAdjacentCell(const int axis, const VectorDi &face, const int ord) { return face - VectorDi::Unit(axis) * (ord & 1 ^ 1); }

	static constexpr int cellEdgeAxis(const int ord) { return ord >> (Dim - 1); }
	static VectorDi edgeAdjacentNode(const int axis, const VectorDi &edge, const int ord) { return edge + VectorDi::Unit(axis) * (ord & 1); }

	static VectorDi cellEdge(const VectorDi &cell, const int ord)
	{
		if constexpr (Dim == 2) return cell + VectorDi::Unit(cellEdgeAxis(ord) ^ 1) * (ord & 1);
		else return cell + VectorDi::Unit((cellEdgeAxis(ord) + 1) % 3) * (ord & 1) + VectorDi::Unit((cellEdgeAxis(ord) + 2) % 3) * (ord >> 1 & 1);
	}
};

}
