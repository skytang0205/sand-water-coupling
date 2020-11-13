#pragma once

#include "Structures/Array.h"
#include "Structures/Field.h"
#include "Structures/Grid.h"

namespace PhysX {

template <int Dim>
class GridBasedScalarField : public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Grid<Dim> *_grid;

public:

	GridBasedScalarField(const Grid<Dim> *const grid) : _grid(grid) { }

	virtual ~GridBasedScalarField() = default;

	const Grid<Dim> &grid() const { return *_grid; }

	virtual real operator()(const VectorDr &pos) const = 0;
	virtual VectorDr gradient(const VectorDr &pos) const = 0;

	virtual void read(std::istream &in) = 0;
	virtual void write(std::ostream &out) const = 0;
};

template <int Dim>
class CellCenteredScalarField final : public GridBasedScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

	using GridBasedScalarField<Dim>::_grid;

protected:

	Array<Dim> _data;

public:

	CellCenteredScalarField(const Grid<Dim> *const grid, const real value = 0) : GridBasedScalarField<Dim>(grid), _data(grid->cellSize(), value) { }

	virtual ~CellCenteredScalarField() = default;

	const real &operator[](const VectorDi &coord) const { return _data[coord]; }
	real &operator[](const VectorDi &coord) { return _data[coord]; }

	virtual real operator()(const VectorDr &pos) const override;

	VectorDr gradientAtCellCenter(const VectorDi &cell) const;
	virtual VectorDr gradient(const VectorDr &pos) const override;

	virtual void read(std::istream &in) { _data.read(in); }
	virtual void write(std::ostream &out) const { _data.write(out); }
};

}
