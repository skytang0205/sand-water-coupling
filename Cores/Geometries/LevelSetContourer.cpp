#include "LevelSetContourer.h"

#include "Structures/StaggeredGrid.h"

#include <iostream>

#include <cstdlib>

namespace PhysX {

template <>
MarchingCubesContourer<2>::MarchingCubesContourer(const Grid<2> *const nodeGrid) :
	_nodeGrid(nodeGrid),
	_cellGrid(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Ones(), _nodeGrid->dataOrigin() + VectorDr::Ones() * _nodeGrid->spacing() / 2),
	_edgeGrids {
		Grid<2>(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Unit(0), _nodeGrid->dataOrigin() + VectorDr::Unit(0) * _nodeGrid->spacing() / 2),
		Grid<2>(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Unit(1), _nodeGrid->dataOrigin() + VectorDr::Unit(1) * _nodeGrid->spacing() / 2)
	},
	_edgeMark {
		GridBasedData<2, int>(&_edgeGrids[0], -1),
		GridBasedData<2, int>(&_edgeGrids[1], -1)
	}
{ }

template <>
MarchingCubesContourer<3>::MarchingCubesContourer(const Grid<3> *const nodeGrid) :
	_nodeGrid(nodeGrid),
	_cellGrid(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Ones(), _nodeGrid->dataOrigin() + VectorDr::Ones() * _nodeGrid->spacing() / 2),
	_edgeGrids {
		Grid<3>(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Unit(0), _nodeGrid->dataOrigin() + VectorDr::Unit(0) * _nodeGrid->spacing() / 2),
		Grid<3>(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Unit(1), _nodeGrid->dataOrigin() + VectorDr::Unit(1) * _nodeGrid->spacing() / 2),
		Grid<3>(_nodeGrid->spacing(), _nodeGrid->dataSize() - VectorDi::Unit(2), _nodeGrid->dataOrigin() + VectorDr::Unit(2) * _nodeGrid->spacing() / 2)
	},
	_edgeMark {
		GridBasedData<3, int>(&_edgeGrids[0], -1),
		GridBasedData<3, int>(&_edgeGrids[1], -1),
		GridBasedData<3, int>(&_edgeGrids[2], -1)
	}
{ }

template <int Dim>
void MarchingCubesContourer<Dim>::contour(const LevelSet<Dim> &levelSet, std::vector<VectorDr> &positions, std::vector<VectorDr> &normals, std::vector<uint> &indices)
{
	const auto &sdf = levelSet.signedDistanceField();
	if (sdf.grid() != _nodeGrid) {
		std::cerr << "Error: [MarchingCubesContourer] encountered incompatible grid." << std::endl;
		std::exit(-1);
	}

	positions.clear();
	normals.clear();
	indices.clear();

	_cellGrid.forEach([&](const VectorDi &cell) {
		const uint cellType = getCellType(sdf, cell);
		uint edgeState;
		if constexpr (Dim == 2) edgeState = _kEdgeStateTable2[cellType];
		else edgeState = _kEdgeStateTable3[cellType];
		for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellEdges(); i++) {
			if (edgeState >> i & 1) {
				const int axis = StaggeredGrid<Dim>::cellEdgeAxis(i);
				const VectorDi edge = StaggeredGrid<Dim>::cellEdge(cell, i);
				if (_edgeMark[axis][edge] < 0) {
					const VectorDi node0 = StaggeredGrid<Dim>::edgeAdjacentNode(axis, edge, 0);
					const VectorDi node1 = StaggeredGrid<Dim>::edgeAdjacentNode(axis, edge, 1);
					const real theta = Surface<Dim>::theta(sdf[node0], sdf[node1]);
					const VectorDr pos = (1 - theta) * sdf.position(node0) + theta * sdf.position(node1);
					_edgeMark[axis][edge] = int(positions.size());
					positions.push_back(pos);
				}
			}
		}
		const int *edgeOrds;
		if constexpr (Dim == 2) edgeOrds = _kEdgeOrdsTable2[cellType];
		else edgeOrds = _kEdgeOrdsTable3[cellType];
		for (int i = 0; edgeOrds[i] != -1; i++) {
			const int axis = StaggeredGrid<Dim>::cellEdgeAxis(edgeOrds[i]);
			const VectorDi edge = StaggeredGrid<Dim>::cellEdge(cell, edgeOrds[i]);
			indices.push_back(uint(_edgeMark[axis][edge]));
		}
	});

	normals.resize(positions.size());
	for (size_t i = 0; i < positions.size(); i++)
		normals[i] = levelSet.closestNormal(positions[i]);
}

template <int Dim>
uint MarchingCubesContourer<Dim>::getCellType(const GridBasedScalarField<Dim> &sdf, const VectorDi &cell) const
{
	uint type = 0;
	for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellNodes(); i++) {
		const VectorDi node = StaggeredGrid<Dim>::cellNode(cell, i);
		if (Surface<Dim>::isInside(sdf[node])) type |= 1 << i;
	}
	// Handle diagonal cases for 2D.
	if constexpr (Dim == 2) {
		if (type == 6 || type == 9) {
			const real phi0 = sdf[StaggeredGrid<Dim>::cellNode(cell, 0)];
			const real phi1 = sdf[StaggeredGrid<Dim>::cellNode(cell, 1)];
			const real phi2 = sdf[StaggeredGrid<Dim>::cellNode(cell, 2)];
			const real phi3 = sdf[StaggeredGrid<Dim>::cellNode(cell, 3)];
			const real centerPhi = (phi0 + phi1 + phi2 + phi3) * real(0.25);
			if (Surface<Dim>::isInside(centerPhi)) type = 16 + (type < 8);
		}
	}
	return type;
}

template class MarchingCubesContourer<2>;
template class MarchingCubesContourer<3>;

}
