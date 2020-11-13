#pragma once

#include "Structures/Array.h"
#include "Structures/Field.h"
#include "Structures/Grid.h"

namespace PhysX {

template <int Dim>
class GridBasedVectorField : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Grid<Dim> *_grid;

public:

	GridBasedVectorField(const Grid<Dim> *const grid) : _grid(grid) { }

	virtual ~GridBasedVectorField() = default;

	const Grid<Dim> &grid() const { return *_grid; }

	virtual VectorDr operator()(const VectorDr &pos) const = 0;
	virtual real divergence(const VectorDr &pos) const = 0;

	virtual void read(std::istream &in) = 0;
	virtual void write(std::ostream &out) const = 0;
};

template <int Dim>
class FaceCenteredVectorField final : public GridBasedVectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

	using GridBasedVectorField<Dim>::_grid;

protected:

	std::array<Array<Dim>, Dim> _data;

public:

	FaceCenteredVectorField(const Grid<Dim> *const grid, const VectorDr &value = VectorDr::Zero()) : GridBasedVectorField<Dim>(grid)
	{
		for (int axis = 0; axis < Dim; axis++) {
			_data[axis].resize(_grid->faceSize(axis), value[axis]);
		}
	}

	virtual ~FaceCenteredVectorField() = default;

	const Array<Dim> &operator[](const int axis) const { return _data[axis]; }
	Array<Dim> &operator[](const int axis) { return _data[axis]; }

	real operator()(const int axis, const VectorDr &pos) const;
	virtual VectorDr operator()(const VectorDr &pos) const override final;

	real divergenceAtCellCenter(const VectorDi &cell) const;
	virtual real divergence(const VectorDr &pos) const override final;

	virtual void read(std::istream &in) { for (int axis = 0; axis < Dim; axis++) _data[axis].read(in); }
	virtual void write(std::ostream &out) const { for (int axis = 0; axis < Dim; axis++) _data[axis].write(out); }
};

}
