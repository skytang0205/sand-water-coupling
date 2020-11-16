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

	virtual VectorDr closestPoint(const VectorDr &pos) const = 0;
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const = 0;
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const = 0;

	real fractionInside(const VectorDr &pos0, const VectorDr &pos1) const
	{
		const real phi0 = signedDistance(pos0);
		const real phi1 = signedDistance(pos1);
		if (phi0 < 0 && phi1 < 0) return 1;
		else if (phi0 < 0 && phi1 >= 0) return phi0 / (phi0 - phi1);
		else if (phi0 >= 0 && phi1 < 0) return phi1 / (phi1 - phi0);
		else return 0;
	}
};

}
