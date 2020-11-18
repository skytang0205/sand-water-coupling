#pragma once

#include "Geometries/GridBasedImplicitSurface.h"

#include <queue>
#include <vector>

namespace PhysX {

template <int Dim>
class LevelSetReinitializer
{
	DECLARE_DIM_TYPES(Dim)

protected:

	int _reinitMaxIters;

public:

	LevelSetReinitializer(const int reinitMaxIters) : _reinitMaxIters(reinitMaxIters) { }

	LevelSetReinitializer(const LevelSetReinitializer &rhs) = delete;
	LevelSetReinitializer &operator=(const LevelSetReinitializer &rhs) = delete;
	virtual ~LevelSetReinitializer() = default;

	virtual void reinitialize(GridBasedImplicitSurface<Dim> &levelSet) = 0;
};

template <int Dim>
class FastMarchingReinitializer final : public LevelSetReinitializer<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using LevelSetReinitializer<Dim>::_reinitMaxIters;

	GridBasedData<Dim> _tent; // tentative signed distance
	GridBasedData<Dim, uchar> _visited;

	std::vector<int> _intfCellIndices;
	std::priority_queue<std::pair<real, int>> _heap;

public:

	FastMarchingReinitializer(const Grid<Dim> *const grid, const int reinitMaxIters);

	virtual void reinitialize(GridBasedImplicitSurface<Dim> &levelSet) override;

protected:

	void initInterface(const GridBasedScalarField<Dim> &phi);
	void performFastMarching();

	void updateNeighborCells(const VectorDi &cell);
	real solveEikonalEquation(const VectorDi &cell) const;

	static real solveQuadratic(const real p0, const real dx) { return p0 + dx; }
	static real solveQuadratic(real p0, real p1, const real dx);
	static real solveQuadratic(real p0, real p1, real p2, const real dx);
};

}
