#include "LevelSetReinitializer.h"

#include <algorithm>
#include <limits>

#include <cmath>

namespace PhysX {

template <int Dim>
FastMarchingReinitializer<Dim>::FastMarchingReinitializer(const Grid<Dim> *const grid, const int reinitMaxIters) :
	LevelSetReinitializer<Dim>(reinitMaxIters),
	_tent(grid),
	_visited(grid)
{ }

template <int Dim>
void FastMarchingReinitializer<Dim>::reinitialize(GridBasedImplicitSurface<Dim>&levelSet)
{
	auto &phi = levelSet.signedDistanceField();
	const real bandWidth = _reinitMaxIters * phi.spacing();
	_tent.setConstant(bandWidth > 0 ? bandWidth : std::numeric_limits<real>::infinity());
	_visited.setZero();
	initInterface(phi);
	performFastMarching();
	phi.parallelForEach([&](const VectorDi &cell) {
		phi[cell] = (phi[cell] > 0 ? 1 : -1) * _tent[cell];
	});
}

template <int Dim>
void FastMarchingReinitializer<Dim>::initInterface(const GridBasedScalarField<Dim> &phi)
{
	const auto isInterface = [](const real phi0, const real phi1)->bool {
		return (phi0 <= 0 && phi1 > 0) || (phi0 > 0 && phi1 <= 0);
	};

	_intfCellIndices.clear();
	phi.forEach([&](const VectorDi &cell) {
		VectorDr tempPhi = VectorDr::Ones() * std::numeric_limits<real>::infinity();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			if (phi.isValid(nbCell) && isInterface(phi[cell], phi[nbCell])) {
				const int axis = i >> 1;
				tempPhi[axis] = std::min(tempPhi[axis], phi[cell] / (phi[cell] - phi[nbCell]));
			}
		}
		if (tempPhi.array().isFinite().any()) {
			_tent[cell] = real(1) / tempPhi.cwiseInverse().norm();
			_visited[cell] = true;
			_intfCellIndices.push_back(int(phi.index(cell)));
		}
	});
}

template <int Dim>
void FastMarchingReinitializer<Dim>::performFastMarching()
{
	for (const auto cellIdx : _intfCellIndices)
		updateNeighborCells(_tent.coordinate(cellIdx));
	while (!_heap.empty()) {
		const real val = _heap.top().first;
		const VectorDi cell = _tent.coordinate(_heap.top().second);
		_heap.pop();
		if (_tent[cell] != val) continue;
		_visited[cell] = true;
		updateNeighborCells(cell);
	}
}

template <int Dim>
void FastMarchingReinitializer<Dim>::updateNeighborCells(const VectorDi &cell)
{
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
		if (!_tent.isValid(nbCell)) continue;
		if (const real temp = solveEikonalEquation(nbCell); temp < _tent[nbCell]) {
			_tent[nbCell] = temp;
			_heap.push(std::make_pair(temp, int(_tent.index(nbCell))));
		}
	}
}

template <int Dim>
real FastMarchingReinitializer<Dim>::solveEikonalEquation(const VectorDi &cell) const
{
	VectorDr tempPhi = VectorDr::Ones() * std::numeric_limits<real>::infinity();
	VectorDi mark = VectorDi::Zero();
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
		if (_tent.isValid(nbCell) && _visited[nbCell]) {
			const int axis = i >> 1;
			tempPhi[axis] = std::min(tempPhi[axis], _tent[nbCell]);
		}
	}
	real newPhi;
	if constexpr (Dim == 2) newPhi = solveQuadratic(tempPhi.x(), tempPhi.y(), _tent.spacing());
	else newPhi = solveQuadratic(tempPhi.x(), tempPhi.y(), tempPhi.z(), _tent.spacing());
	if (!std::isfinite(newPhi)) {
		std::cerr << "Error: [FastMarchingReinitializer] failed to solve Eikonal equation." << std::endl;
		std::exit(-1);
	}
	return newPhi;
}

template <int Dim>
real PhysX::FastMarchingReinitializer<Dim>::solveQuadratic(real p0, real p1, const real dx)
{
	if (p0 > p1) std::swap(p0, p1);
	if (std::isinf(p1) && p1 - p0 > dx) return solveQuadratic(p0, dx);
	else return ((p0 + p1) + std::sqrt(2 * dx * dx - (p0 - p1) * (p0 - p1))) * real(0.5);
}

template <int Dim>
real PhysX::FastMarchingReinitializer<Dim>::solveQuadratic(real p0, real p1, real p2, const real dx)
{
	if (p0 > p1) std::swap(p0, p1);
	if (p1 > p2) std::swap(p1, p2);
	if (std::isinf(p2) && (p2 - p0) * (p2 - p0) + (p2 - p1) * (p2 - p1) > dx * dx) return solveQuadratic(p0, p1, dx);
	else return (p0 + p1 + p2 + std::sqrt((p0 + p1 + p2) * (p0 + p1 + p2) - 3 * (p0 * p0 + p1 * p1 + p2 * p2 - dx * dx))) / 3;
}

template class LevelSetReinitializer<2>;
template class LevelSetReinitializer<3>;

template class FastMarchingReinitializer<2>;
template class FastMarchingReinitializer<3>;

}
