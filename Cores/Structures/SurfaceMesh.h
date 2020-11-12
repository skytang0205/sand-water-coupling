#pragma once

#include "Structures/Geometry.h"

#include <vector>

namespace PhysX {

template <int Dim>
class SurfaceMesh final : public Geometry<Dim>
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

public:

	std::vector<VectorDr> positions;
	std::vector<VectorDr> normals;
	std::vector<uint> indices;

	SurfaceMesh() = default;
	SurfaceMesh(const SurfaceMesh &rhs) = default;
	SurfaceMesh &operator=(const SurfaceMesh &rhs) = default;
	virtual ~SurfaceMesh() = default;
};

}
