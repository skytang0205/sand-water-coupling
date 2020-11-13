#include "GridBasedVectorField.h"

namespace PhysX {

template<int Dim>
real FaceCenteredVectorField<Dim>::operator()(const int axis, const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getFaceLerp(axis, pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += _data[axis][coords[i]] * weights[i];
	return val;
}

template <int Dim>
Vector<real, Dim> FaceCenteredVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec;
	for (int axis = 0; axis < Dim; axis++)
		vec[axis] = operator()(axis, pos);
	return vec;
}

template <int Dim>
real FaceCenteredVectorField<Dim>::divergenceAtCellCenter(const VectorDi &cell) const
{
	real val = 0;
	for (int i = 0; i < Dim; i++) {
		val += _data[i][cell + VectorDi::Unit(i)] - _data[i][cell];
	}
	return val / (2 * _grid->spacing());
}

template <int Dim>
real FaceCenteredVectorField<Dim>::divergence(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getCellLerp(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += divergenceAtCellCenter(coords[i]) * weights[i];
	return val;
}

template class GridBasedVectorField<2>;
template class GridBasedVectorField<3>;

template class FaceCenteredVectorField<2>;
template class FaceCenteredVectorField<3>;
}
