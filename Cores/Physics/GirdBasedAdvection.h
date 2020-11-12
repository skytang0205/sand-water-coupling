#pragma once

#include "Structures/ScalarGridField.h"
#include "Structures/VectorGridField.h"

namespace PhysX {

template <int Dim>
class GridBasedAdvection
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

public:

	GridBasedAdvection() = default;
	GridBasedAdvection(const GridBasedAdvection &rhs) = default;
	GridBasedAdvection &operator=(const GridBasedAdvection &rhs) = default;
	virtual ~GridBasedAdvection() = default;

	void advect(VectorGridField<Dim, FaceCentered> &field, const VectorField<Dim> &flow, const real dt);
};

}
