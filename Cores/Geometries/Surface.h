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
};

}
