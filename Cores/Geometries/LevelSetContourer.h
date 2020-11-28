#pragma once

#include "Geometries/LevelSet.h"

namespace PhysX {

template <int Dim>
class LevelSetContourer
{
	DECLARE_DIM_TYPES(Dim)

public:

	LevelSetContourer() = default;
	LevelSetContourer(const LevelSetContourer &rhs) = delete;
	LevelSetContourer &operator=(const LevelSetContourer &rhs) = delete;
	virtual ~LevelSetContourer() = default;

	virtual void contour(
		const LevelSet<Dim> &levelSet,
		std::vector<VectorDr> &positions,
		std::vector<VectorDr> &normals,
		std::vector<uint> &indices) = 0;
};

template <int Dim>
class MarchingCubesContourer final : public LevelSetContourer<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

#include "MarchingCubesTables.inc"

	const Grid<Dim> *const _nodeGrid;
	const Grid<Dim> _cellGrid;
	const std::array<Grid<Dim>, Dim> _edgeGrids;
	std::array<GridBasedData<Dim, int>, Dim> _edgeMark;

public:

	MarchingCubesContourer(const Grid<Dim> *const nodeGrid);
	MarchingCubesContourer(const MarchingCubesContourer &rhs) = delete;
	MarchingCubesContourer &operator=(const MarchingCubesContourer &rhs) = delete;
	virtual ~MarchingCubesContourer() = default;

	virtual void contour(
		const LevelSet<Dim> &levelSet,
		std::vector<VectorDr> &positions,
		std::vector<VectorDr> &normals,
		std::vector<uint> &indices) override;

protected:

	uint getCellType(const GridBasedScalarField<Dim> &sdf, const VectorDi &cell) const;
};

}
