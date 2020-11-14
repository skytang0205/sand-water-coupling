#include "StaggeredGridBasedVectorField.h"

namespace PhysX {

template <>
StaggeredGridBasedVectorField<2>::StaggeredGridBasedVectorField(const StaggeredGrid<2> *grid, const VectorDr &value) :
	_grid(grid),
	_components {
		GridBasedScalarField<2>(_grid->faceGrid(0), value[0]),
		GridBasedScalarField<2>(_grid->faceGrid(1), value[1])
}
{ }

template <>
StaggeredGridBasedVectorField<3>::StaggeredGridBasedVectorField(const StaggeredGrid<3> *grid, const VectorDr &value) :
	_grid(grid),
	_components {
		GridBasedScalarField<3>(_grid->faceGrid(0), value[0]),
		GridBasedScalarField<3>(_grid->faceGrid(1), value[1]),
		GridBasedScalarField<3>(_grid->faceGrid(2), value[2])
}
{ }

template <int Dim>
Vector<real, Dim> StaggeredGridBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec;
	for (int axis = 0; axis < Dim; axis++)
		vec[axis] = _components[axis](pos);
	return vec;
}

template <int Dim>
real StaggeredGridBasedVectorField<Dim>::divergenceAtCellCenter(const VectorDi &cell) const
{
	real val = 0;
	for (int axis = 0; axis < Dim; axis++) {
		val += _components[axis][cell + VectorDi::Unit(axis)] - _components[axis][cell];
	}
	return val / (2 * _grid->spacing());
}

template <int Dim>
real StaggeredGridBasedVectorField<Dim>::divergence(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->cellGrid()->getLerpCoordsAndWeights(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += divergenceAtCellCenter(coords[i]) * weights[i];
	return val;
}

template class StaggeredGridBasedVectorField<2>;
template class StaggeredGridBasedVectorField<3>;

}
