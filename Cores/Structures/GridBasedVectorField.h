#pragma once

#include "Structures/Field.h"
#include "Structures/GridBasedData.h"

namespace PhysX {

template <int Dim>
class GridBasedVectorField final : public VectorField<Dim>, public GridBasedVectorData<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using GridBasedVectorData<Dim>::_grid;

public:

	using GridBasedVectorData<Dim>::at;

	GridBasedVectorField(const Grid<Dim> *const grid, const VectorDr &value = VectorDr::Zero()) : GridBasedVectorData<Dim>(grid, value) { }

	GridBasedVectorField() = default;
	virtual ~GridBasedVectorField() = default;

	virtual VectorDr operator()(const VectorDr &pos) const override;
	real divergenceAtDataPoint(const VectorDi &coord) const;
	virtual real divergence(const VectorDr &pos) const override;
};

}
