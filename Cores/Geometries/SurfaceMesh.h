#pragma once

#include "Geometries/GridBasedImplicitSurface.h"

#include <vector>

namespace PhysX {

template <int Dim>
class SurfaceMesh final : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	std::vector<VectorDr> positions;
	std::vector<VectorDr> normals;
	std::vector<uint> indices;

	SurfaceMesh() = default;
	SurfaceMesh(const GridBasedImplicitSurface<Dim> &levelSet);
	virtual ~SurfaceMesh() = default;
};

}
