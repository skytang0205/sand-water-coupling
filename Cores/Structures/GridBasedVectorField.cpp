#include "GridBasedVectorField.h"

namespace PhysX {

template <>
GridBasedStaggeredVectorField<2>::GridBasedStaggeredVectorField(const StaggeredGrid<2> *staggeredGrid, const VectorDr &value) :
	_staggeredGrid(staggeredGrid),
	_components {
		GridBasedScalarField<2>(_staggeredGrid->faceGrid(0), value[0]),
		GridBasedScalarField<2>(_staggeredGrid->faceGrid(1), value[1])
	}
{ }

template <>
GridBasedStaggeredVectorField<3>::GridBasedStaggeredVectorField(const StaggeredGrid<3> *staggeredGrid, const VectorDr &value) :
	_staggeredGrid(staggeredGrid),
	_components {
		GridBasedScalarField<3>(_staggeredGrid->faceGrid(0), value[0]),
		GridBasedScalarField<3>(_staggeredGrid->faceGrid(1), value[1]),
		GridBasedScalarField<3>(_staggeredGrid->faceGrid(2), value[2])
	}
{ }

template <int Dim>
Vector<real, Dim> GridBasedStaggeredVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec;
	for (int axis = 0; axis < Dim; axis++)
		vec[axis] = _components[axis](pos);
	return vec;
}

template <int Dim>
real GridBasedStaggeredVectorField<Dim>::divergenceAtCellCenter(const VectorDi &cell) const
{
	real val = 0;
	for (int axis = 0; axis < Dim; axis++) {
		val += _components[axis][cell + VectorDi::Unit(axis)] - _components[axis][cell];
	}
	return val / (2 * _staggeredGrid->spacing());
}

template <int Dim>
real GridBasedStaggeredVectorField<Dim>::divergence(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_staggeredGrid->cellGrid()->getLerpCoordsAndWeights(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += divergenceAtCellCenter(coords[i]) * weights[i];
	return val;
}

template class GridBasedVectorField<2>;
template class GridBasedVectorField<3>;

template class GridBasedStaggeredVectorField<2>;
template class GridBasedStaggeredVectorField<3>;

}
