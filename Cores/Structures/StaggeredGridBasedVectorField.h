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

	const StaggeredGrid<Dim> *_grid = nullptr;
	std::array<GridBasedScalarField<Dim>, Dim> _components;

public:

	StaggeredGridBasedVectorField(const StaggeredGrid<Dim> *const grid, const VectorDr &value = VectorDr::Zero()) { resize(grid, value); }

	StaggeredGridBasedVectorField() = default;
	virtual ~StaggeredGridBasedVectorField() = default;

	void resize(const StaggeredGrid<Dim> *const grid, const VectorDr &value = VectorDr::Zero())
	{
		_grid = grid;
		for (int axis = 0; axis < Dim; axis++)
			_components[axis].resize(_grid->faceGrid(axis), value[axis]);
	}

	bool isBoundary(const int axis, const VectorDi &face) { return _grid->isBoundaryFace(axis, face); }

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
