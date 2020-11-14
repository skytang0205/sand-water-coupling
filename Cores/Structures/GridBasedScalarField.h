#pragma once

#include "Structures/Field.h"
#include "Structures/GridBasedData.h"

namespace PhysX {

template <int Dim>
class GridBasedScalarField final : public ScalarField<Dim>, public GridBasedData<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using GridBasedData<Dim>::_grid;

public:

	using GridBasedData<Dim>::operator[];

	GridBasedScalarField(const Grid<Dim> *const grid, const real value = 0) : GridBasedData<Dim>(grid, value) { }

	virtual ~GridBasedScalarField() = default;

	virtual real operator()(const VectorDr &pos) const override final;
	VectorDr gradientAtDataPoint(const VectorDi &coord) const;
	virtual VectorDr gradient(const VectorDr &pos) const override final;
};

}
