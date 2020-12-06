#pragma once

#include "Geometries/LevelSet.h"

#include <queue>
#include <vector>

namespace PhysX {

template <int Dim>
class LevelSetReinitializer
{
	DECLARE_DIM_TYPES(Dim)

public:

	LevelSetReinitializer() = default;
	LevelSetReinitializer(const LevelSetReinitializer &rhs) = delete;
	LevelSetReinitializer &operator=(const LevelSetReinitializer &rhs) = delete;
	virtual ~LevelSetReinitializer() = default;

	virtual void reinitialize(LevelSet<Dim> &levelSet, const int maxSteps = -1) = 0;
};

template <int Dim>
class FastMarchingReinitializer final : public LevelSetReinitializer<Dim>
{
	DECLARE_DIM_TYPES(Dim)

	using HeapElement = std::pair<real, int>;

protected:

	GridBasedScalarData<Dim> _tent; // tentative signed distance
	GridBasedData<Dim, uchar> _visited;

	std::vector<int> _intfIndices;
	std::priority_queue<HeapElement, std::vector<HeapElement>, std::greater<HeapElement>> _heap;

public:

	FastMarchingReinitializer(const Grid<Dim> *const grid);
	FastMarchingReinitializer(const FastMarchingReinitializer &rhs) = delete;
	FastMarchingReinitializer &operator=(const FastMarchingReinitializer &rhs) = delete;
	virtual ~FastMarchingReinitializer() = default;

	virtual void reinitialize(LevelSet<Dim> &levelSet, const int maxSteps = -1) override;

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
