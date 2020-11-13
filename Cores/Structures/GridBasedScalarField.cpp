#include "GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
real CellCenteredScalarField<Dim>::operator()(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getCellLerp(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += _data[coords[i]] * weights[i];
	return val;
}

template <int Dim>
Vector<real, Dim> CellCenteredScalarField<Dim>::gradientAtCellCenter(const VectorDi &cell) const
{
	VectorDr vec;
	for (int i = 0; i < Dim; i++) {
		vec[i] = _data[cell[i] < _grid->cellSize()[i] - 1 ? cell + VectorDi::Unit(i) : cell]
			- _data[cell[i] > 0 ? cell - VectorDi::Unit(i) : cell];
	}
	return vec / (2 * _grid->spacing());
}

template <int Dim>
Vector<real, Dim> CellCenteredScalarField<Dim>::gradient(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getCellLerp(pos, coords, weights);

	VectorDr vec = VectorDr::Zero();
	for (int i = 0; i < (1 << Dim); i++)
		vec += gradientAtCellCenter(coords[i]) * weights[i];
	return vec;
}

template class GridBasedScalarField<2>;
template class GridBasedScalarField<3>;

template class CellCenteredScalarField<2>;
template class CellCenteredScalarField<3>;

}
