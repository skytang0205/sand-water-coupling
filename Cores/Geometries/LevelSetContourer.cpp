#include "LevelSetContourer.h"

namespace PhysX {

void MarchingCubesContourer<2>::contour(
	const GridBasedImplicitSurface<2> &levelSet,
	std::vector<Vector2r> &positions,
	std::vector<Vector2r> &normals,
	std::vector<uint> &indicess,
	const real isoValue)
{
}

void MarchingCubesContourer<3>::contour(
	const GridBasedImplicitSurface<3> &levelSet,
	std::vector<Vector3r> &positions,
	std::vector<Vector3r> &normals,
	std::vector<uint> &indicess,
	const real isoValue)
{
}

}
