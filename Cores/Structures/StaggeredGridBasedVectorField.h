#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/IO.h"

namespace PhysX {

template <int Dim>
class StaggeredGridBasedVectorField final : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> *_grid;
	std::array<GridBasedScalarField<Dim>, Dim> _components;

public:

	StaggeredGridBasedVectorField(const StaggeredGrid<Dim> *grid, const VectorDr &value = VectorDr::Zero());

	virtual ~StaggeredGridBasedVectorField() = default;

	const GridBasedScalarField<Dim> &operator[](const int axis) const { return _components[axis]; }
	GridBasedScalarField<Dim> &operator[](const int axis) { return _components[axis]; }

	virtual VectorDr operator()(const VectorDr &pos) const override;

	real divergenceAtCellCenter(const VectorDi &cell) const;
	virtual real divergence(const VectorDr &pos) const override;

	void forEach(const std::function<void(const int, const VectorDi &)> &func) const { _grid->forEachFace(func); }
	void parallelForEach(const std::function<void(const int, const VectorDi &)> &func) const { _grid->parallelForEachFace(func); }

	void read(std::istream &in) { for (int axis = 0; axis < Dim; axis++) _components[axis].read(in); }
	void write(std::ostream &out) const { for (int axis = 0; axis < Dim; axis++) _components[axis].write(out); }
};

}
