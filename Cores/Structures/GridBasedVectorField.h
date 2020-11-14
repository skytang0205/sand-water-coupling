#pragma once

#include "Structures/Field.h"
#include "Structures/GridBasedScalarField.h"
#include "Structures/StaggeredGrid.h"

namespace PhysX {

template <int Dim>
class GridBasedVectorField : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	GridBasedVectorField() = default;
	virtual ~GridBasedVectorField() = default;

	virtual VectorDr operator()(const VectorDr &pos) const = 0;
	virtual real divergence(const VectorDr &pos) const = 0;

	virtual void read(std::istream &in) = 0;
	virtual void write(std::ostream &out) const = 0;
};

template <int Dim>
class GridBasedStaggeredVectorField final : public GridBasedVectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> *_staggeredGrid;
	std::array<GridBasedScalarField<Dim>, Dim> _components;

public:

	GridBasedStaggeredVectorField(const StaggeredGrid<Dim> *staggeredGrid, const VectorDr &value = VectorDr::Zero());

	virtual ~GridBasedStaggeredVectorField() = default;

	const GridBasedScalarField<Dim> &operator[](const int axis) const { return _components[axis]; }
	GridBasedScalarField<Dim> &operator[](const int axis) { return _components[axis]; }

	virtual VectorDr operator()(const VectorDr &pos) const override final;

	real divergenceAtCellCenter(const VectorDi &cell) const;
	virtual real divergence(const VectorDr &pos) const override final;

	void forEach(const std::function<void(const int, const VectorDi &)> &func) const { _staggeredGrid->forEachFace(func); }
	void parallelForEach(const std::function<void(const int, const VectorDi &)> &func) const { _staggeredGrid->parallelForEachFace(func); }

	virtual void read(std::istream &in) override { for (int axis = 0; axis < Dim; axis++) _components[axis].read(in); }
	virtual void write(std::ostream &out) const override { for (int axis = 0; axis < Dim; axis++) _components[axis].write(out); }
};

}
