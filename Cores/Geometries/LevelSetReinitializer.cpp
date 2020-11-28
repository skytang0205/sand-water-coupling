#include "LevelSetReinitializer.h"

#include <algorithm>
#include <limits>

#include <cmath>

namespace PhysX {

template <int Dim>
FastMarchingReinitializer<Dim>::FastMarchingReinitializer(const Grid<Dim> *const grid) :
	_tent(grid),
	_visited(grid)
{ }

template <int Dim>
void FastMarchingReinitializer<Dim>::reinitialize(GridBasedImplicitSurface<Dim> &levelSet, const int maxSteps)
{
	auto &phi = levelSet.signedDistanceField();
	const real bandWidth = maxSteps * phi.spacing();
	_tent.setConstant(bandWidth > 0 ? bandWidth : std::numeric_limits<real>::infinity());
	_visited.setZero();
	initInterface(phi);
	performFastMarching();
	phi.parallelForEach([&](const VectorDi &coord) {
		phi[coord] = Surface<Dim>::sign(phi[coord]) * _tent[coord];
	});
}

template <int Dim>
void FastMarchingReinitializer<Dim>::initInterface(const GridBasedScalarField<Dim> &phi)
{
	_intfIndices.clear();
	phi.forEach([&](const VectorDi &coord) {
		VectorDr tempPhi = VectorDr::Ones() * std::numeric_limits<real>::infinity();
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
			if (phi.isValid(nbCoord) && Surface<Dim>::isInterface(phi[coord], phi[nbCoord])) {
				const int axis = Grid<Dim>::neighborAxis(i);
				tempPhi[axis] = std::min(tempPhi[axis], Surface<Dim>::theta(phi[coord], phi[nbCoord]) * phi.spacing());
			}
		}
		if (tempPhi.array().isFinite().any()) {
			_tent[coord] = real(1) / tempPhi.cwiseInverse().norm();
			_visited[coord] = true;
			_intfIndices.push_back(int(phi.index(coord)));
		}
	});
}

template <int Dim>
void FastMarchingReinitializer<Dim>::performFastMarching()
{
	for (const auto index : _intfIndices)
		updateNeighbors(_tent.coordinate(index));
	while (!_heap.empty()) {
		const real val = _heap.top().first;
		const VectorDi coord = _tent.coordinate(_heap.top().second);
		_heap.pop();
		if (_tent[coord] != val) continue;
		_visited[coord] = true;
		updateNeighbors(coord);
	}
}

template <int Dim>
void FastMarchingReinitializer<Dim>::updateNeighbors(const VectorDi &coord)
{
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
		if (!_tent.isValid(nbCoord) || _visited[nbCoord]) continue;
		if (const real temp = solveEikonalEquation(nbCoord); temp < _tent[nbCoord]) {
			_tent[nbCoord] = temp;
			_heap.push(HeapElement(temp, int(_tent.index(nbCoord))));
		}
	}
}

template <int Dim>
real FastMarchingReinitializer<Dim>::solveEikonalEquation(const VectorDi &coord) const
{
	VectorDr tempPhi = VectorDr::Ones() * std::numeric_limits<real>::infinity();
	VectorDi mark = VectorDi::Zero();
	for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
		const VectorDi nbCoord = Grid<Dim>::neighbor(coord, i);
		if (_tent.isValid(nbCoord) && _visited[nbCoord]) {
			const int axis = Grid<Dim>::neighborAxis(i);
			tempPhi[axis] = std::min(tempPhi[axis], _tent[nbCoord]);
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
	if (std::isinf(p1) || p1 - p0 > dx) return solveQuadratic(p0, dx);
	else return ((p0 + p1) + std::sqrt(2 * dx * dx - (p0 - p1) * (p0 - p1))) * real(0.5);
}

template <int Dim>
real PhysX::FastMarchingReinitializer<Dim>::solveQuadratic(real p0, real p1, real p2, const real dx)
{
	if (p0 > p1) std::swap(p0, p1);
	if (p1 > p2) std::swap(p1, p2);
	if (std::isinf(p2) || (p2 - p0) * (p2 - p0) + (p2 - p1) * (p2 - p1) > dx * dx) return solveQuadratic(p0, p1, dx);
	else return (p0 + p1 + p2 + std::sqrt((p0 + p1 + p2) * (p0 + p1 + p2) - 3 * (p0 * p0 + p1 * p1 + p2 * p2 - dx * dx))) / 3;
}

template class LevelSetReinitializer<2>;
template class LevelSetReinitializer<3>;

template class FastMarchingReinitializer<2>;
template class FastMarchingReinitializer<3>;

}
