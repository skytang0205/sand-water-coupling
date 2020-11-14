#include "StaggeredGridBasedVectorField.h"

namespace PhysX {

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
