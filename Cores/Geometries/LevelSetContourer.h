#pragma once

#include "Geometries/GridBasedImplicitSurface.h"

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
		const GridBasedImplicitSurface<Dim> &levelSet,
		std::vector<VectorDr> &positions,
		std::vector<VectorDr> &normals,
		std::vector<uint> &indices,
		const real isoValue = 0) = 0;
};

template <int Dim> class MarchingCubesContourer;

template <>
class MarchingCubesContourer<2> : public LevelSetContourer<2>
{
public:

	MarchingCubesContourer() = default;
	MarchingCubesContourer(const MarchingCubesContourer &rhs) = delete;
	MarchingCubesContourer &operator=(const MarchingCubesContourer &rhs) = delete;
	virtual ~MarchingCubesContourer() = default;

	virtual void contour(
		const GridBasedImplicitSurface<2> &levelSet,
		std::vector<Vector2r> &positions,
		std::vector<Vector2r> &normals,
		std::vector<uint> &indicess,
		const real isoValue = 0) override;
};

template <>
class MarchingCubesContourer<3> : public LevelSetContourer<3>
{
public:

	MarchingCubesContourer() = default;
	MarchingCubesContourer(const MarchingCubesContourer &rhs) = delete;
	MarchingCubesContourer &operator=(const MarchingCubesContourer &rhs) = delete;
	virtual ~MarchingCubesContourer() = default;

	virtual void contour(
		const GridBasedImplicitSurface<3> &levelSet,
		std::vector<Vector3r> &positions,
		std::vector<Vector3r> &normals,
		std::vector<uint> &indicess,
		const real isoValue = 0) override;
};

}
