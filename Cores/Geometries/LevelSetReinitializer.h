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
	using PRI = std::pair<real, int>;

	GridBasedData<Dim> _tent; // tentative signed distance
	GridBasedData<Dim, uchar> _visited;

	std::vector<int> _intfIndices;
	std::priority_queue<PRI, std::vector<PRI>, std::greater<PRI>> _heap;

public:

	FastMarchingReinitializer(const Grid<Dim> *const grid, const int reinitMaxIters);

	virtual void reinitialize(GridBasedImplicitSurface<Dim> &levelSet) override;

protected:

	void initInterface(const GridBasedScalarField<Dim> &phi);
	void performFastMarching();

	void updateNeighbors(const VectorDi &coord);
	real solveEikonalEquation(const VectorDi &coord) const;

	static real solveQuadratic(const real p0, const real dx) { return p0 + dx; }
	static real solveQuadratic(real p0, real p1, const real dx);
	static real solveQuadratic(real p0, real p1, real p2, const real dx);
};

}
