#pragma once

#include "Structures/Geometry.h"

#include <vector>

namespace PhysX {

template <int Dim>
class SurfaceMesh final : public Geometry<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	std::vector<VectorDr> positions;
	std::vector<VectorDr> normals;
	std::vector<uint> indices;

	SurfaceMesh() = default;
	virtual ~SurfaceMesh() = default;
};

}
