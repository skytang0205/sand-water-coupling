#include "SurfaceMesh.h"

#include "Geometries/LevelSetContourer.h"

namespace PhysX {

template <int Dim>
SurfaceMesh<Dim>::SurfaceMesh(const LevelSet<Dim> &levelSet)
{
	std::unique_ptr<LevelSetContourer<Dim>> contourer = std::make_unique<MarchingCubesContourer<Dim>>(levelSet.signedDistanceField().grid());
	contourer->contour(levelSet, positions, normals, indices);
}

template class SurfaceMesh<2>;
template class SurfaceMesh<3>;

}
