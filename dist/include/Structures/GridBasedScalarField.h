#pragma once

#include "Structures/Field.h"
#include "Structures/GridBasedData.h"

namespace PhysX {

template <int Dim>
class GridBasedScalarField final : public ScalarField<Dim>, public GridBasedScalarData<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using GridBasedScalarData<Dim>::_grid;

public:

	using GridBasedScalarData<Dim>::at;

	GridBasedScalarField(const Grid<Dim> *const grid, const real value = 0) : GridBasedScalarData<Dim>(grid, value) { }

	GridBasedScalarField() = default;
	virtual ~GridBasedScalarField() = default;

	virtual real operator()(const VectorDr &pos) const override;
	VectorDr gradientAtDataPoint(const VectorDi &coord) const;
	virtual VectorDr gradient(const VectorDr &pos) const override;
	real laplacianAtDataPoint(const VectorDi &coord) const;
	virtual real laplacian(const VectorDr &pos) const override;
};

}
