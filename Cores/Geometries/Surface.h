#pragma once

#include "Utilities/Types.h"

namespace PhysX {

template <int Dim>
class Surface
{
	DECLARE_DIM_TYPES(Dim)

public:

	Surface() = default;
	virtual ~Surface() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const = 0;
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const = 0;
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const = 0;

	static bool inside(const real phi) { return phi <= 0; }
	static int sgn(const real phi) { return inside(phi) ? -1 : 1; }
	static bool isInterface(const real phi0, const real phi1) { return sgn(phi0) != sgn(phi1); }
	static real theta(const real phi0, const real phi1) { return phi0 / (phi0 - phi1); }
	static real fraction(const real phi0, const real phi1)
	{
		if (inside(phi0) && inside(phi1)) return 1;
		else if (inside(phi0) && !inside(phi1)) return theta(phi0, phi1);
		else if (!inside(phi0) && inside(phi1)) return theta(phi1, phi0);
		else return 0;
	}
};

}
