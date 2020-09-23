#pragma once

#include "Types.h"

template <int Dim>
class Surface
{
private:

	DECLARE_DIM_TYPES(Dim)

public:

	Surface() = default;
	virtual ~Surface() = default;

	virtual VectorDr Closest_Point(const VectorDr &pos) const = 0;
	virtual VectorDr Closest_Normal(const VectorDr &pos) const = 0;

	virtual bool Inside(const VectorDr &pos) const = 0;

	virtual double Distance(const VectorDr &pos) const = 0;
	virtual double Signed_Distance(const VectorDr &pos) const = 0;
};
