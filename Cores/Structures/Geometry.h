#pragma once

#include "Utilities/Types.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class Geometry
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

public:

	Geometry() = default;
	Geometry(const Geometry &rhs) = default;
	Geometry &operator=(const Geometry &rhs) = default;
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
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

public:

	ImplicitGeometry() = default;
	ImplicitGeometry(const ImplicitGeometry &rhs) = default;
	ImplicitGeometry &operator=(const ImplicitGeometry &rhs) = default;
	virtual ~ImplicitGeometry() = default;

	virtual VectorDr closestPoint(const VectorDr &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const { return std::abs(signedDistance(pos)); }
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const { return signedDistance(pos) < 0; }
};

}
