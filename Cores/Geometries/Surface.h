#pragma once

#include "Utilities/Types.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class Surface
{
	DECLARE_DIM_TYPES(Dim)

public:

	Surface() = default;
	virtual ~Surface() = default;

	virtual VectorDr closestPoint(const VectorDr &pos) const = 0;
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const = 0;
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const = 0;
};

template <int Dim>
class ImplicitSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	ImplicitSurface() = default;
	virtual ~ImplicitSurface() = default;

	virtual VectorDr closestPoint(const VectorDr &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const { return std::abs(signedDistance(pos)); }
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const { return signedDistance(pos) < 0; }
};

}
