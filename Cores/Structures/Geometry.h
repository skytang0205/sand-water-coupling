#pragma once

#include "Utilities/Types.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class Geometry
{
	DECLARE_DIM_TYPES(Dim)

public:

	Geometry() = default;
	virtual ~Geometry() = default;

	virtual VectorDr closestPoint(const VectorDr &pos) const = 0;
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const = 0;
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const = 0;
};

template <int Dim>
class ImplicitGeometry : public Geometry<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	ImplicitGeometry() = default;
	virtual ~ImplicitGeometry() = default;

	virtual VectorDr closestPoint(const VectorDr &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const { return std::abs(signedDistance(pos)); }
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const { return signedDistance(pos) < 0; }
};

}
