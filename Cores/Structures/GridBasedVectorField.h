#pragma once

#include "Structures/Field.h"
#include "Structures/GridBasedData.h"

namespace PhysX {

template <int Dim>
class GridBasedVectorField final : public VectorField<Dim>, public GridBasedData<Dim, Vector<real, Dim>>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	using GridBasedData<Dim, VectorDr>::_grid;

public:

	using GridBasedData<Dim, VectorDr>::operator[];

	GridBasedVectorField(const Grid<Dim> *const grid, const VectorDr &value = VectorDr::Zero()) : GridBasedData<Dim, VectorDr>(grid, value) { }

	virtual ~GridBasedVectorField() = default;

	virtual VectorDr operator()(const VectorDr &pos) const override;
	real divergenceAtDataPoint(const VectorDi &coord) const;
	virtual real divergence(const VectorDr &pos) const override;
};

}
