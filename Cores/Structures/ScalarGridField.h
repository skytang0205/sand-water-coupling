#pragma once

#include "Structures/Array.h"
#include "Structures/Field.h"
#include "Structures/Grid.h"

namespace PhysX {

template <int Dim, int DataLayout> class ScalarGridField;

template <int Dim>
class ScalarGridField<Dim, CellCentered> : public Grid<Dim>, public ScalarField<Dim>
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

	const real &operator[](const VectorDi &coord) const { return _data[coord]; }
	real &operator[](const VectorDi &coord) { return _data[coord]; }

	virtual real operator()(const VectorDr &pos) const override final
	{
		VectorDr frac;
		VectorDi cell = this->getLowerCell(pos, frac);
		return _data.lerp(cell, frac);
	}

	virtual VectorDr gradient(const VectorDr &pos) const override final
	{
		const real dx = this->_spacing;
		VectorDr vec;
		for (int i = 0; i < Dim; i++) {
			vec[i] = (operator()(pos + VectorDr::Unit(i) * dx) - operator()(pos - VectorDr::Unit(i) * dx)) / (2 * dx);
		}
		return vec;
	}
};

}
