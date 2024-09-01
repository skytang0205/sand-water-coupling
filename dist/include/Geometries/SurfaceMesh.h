#pragma once

#include "Geometries/LevelSet.h"
#include "Utilities/IO.h"

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
	SurfaceMesh(const LevelSet<Dim> &levelSet);
	virtual ~SurfaceMesh() = default;

	// TODO: implement following functions.
	virtual VectorDr closestPosition(const VectorDr &pos) const override { return VectorDr::Zero(); }
	virtual VectorDr closestNormal(const VectorDr &pos) const override { return VectorDr::Zero(); }
	virtual real distance(const VectorDr &pos) const override { return 0; }
	virtual real signedDistance(const VectorDr &pos) const override { return 0; }
	virtual bool isInside(const VectorDr &pos) const override { return false; }
};

}
