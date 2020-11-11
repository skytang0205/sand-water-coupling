#pragma once

#include "Utilities/Types.h"

namespace PhysX {

template <int Dim>
class ScalarField
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

public:

	ScalarField() = default;
	ScalarField(const ScalarField &rhs) = default;
	ScalarField &operator=(const ScalarField &rhs) = default;
	virtual ~ScalarField() = default;

	virtual real operator()(const VectorDr &pos) = 0;
	virtual VectorDr gradient(const VectorDr &pos) = 0;
};

template <int Dim>
class VectorField
{
	DECLARE_DIM_TYPES(Dim)

public:

	VectorField() = default;
	VectorField(const VectorField &rhs) = default;
	VectorField &operator=(const VectorField &rhs) = default;
	~VectorField() = default;

	virtual VectorDr operator()(const VectorDr &pos) = 0;
	virtual real divergence(const VectorDr &pos) = 0;
};

}
