#pragma once

#include "Structures/Array.h"
#include "Structures/Field.h"
#include "Structures/Grid.h"

namespace PhysX {

template <int Dim, int DataLayout> class ScalarGridField;

template <int Dim>
class ScalarGridField<Dim, CellCentered> : public Grid<Dim>, ScalarField<Dim>
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

protected:

	Array<Dim> _data;

public:

	ScalarGridField(
		const real spacing,
		const VectorDi &resolution,
		const VectorDr &center = VectorDr::Zero(),
		const real value = 0)
		:
		Grid<Dim>(spacing, resolution, center),
		_data(resolution, value)
	{ }

	ScalarGridField(const ScalarGridField &rhs) = default;
	ScalarGridField &operator=(const ScalarGridField &rhs) = default;
	virtual ~ScalarGridField() = default;

	virtual real operator()(const VectorDr &pos) const override
	{
		if constexpr (Dim == 2) {
		}
		else {
		}
	}

	virtual VectorDr gradient(const VectorDr &pos) const override { return pos; }
};

}
