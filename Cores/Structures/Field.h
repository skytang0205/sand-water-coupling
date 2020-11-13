#pragma once

#include "Utilities/Types.h"

namespace PhysX {

template <int Dim>
class ScalarField
{
	DECLARE_DIM_TYPES(Dim)

public:

	ScalarField() = default;
	virtual ~ScalarField() = default;

	virtual real operator()(const VectorDr &pos) const = 0;
	virtual VectorDr gradient(const VectorDr &pos) const = 0;
};

template <int Dim>
class VectorField
{
	DECLARE_DIM_TYPES(Dim)

public:

	VectorField() = default;
	~VectorField() = default;

	virtual VectorDr operator()(const VectorDr &pos) const = 0;
	virtual real divergence(const VectorDr &pos) const = 0;
};

}
