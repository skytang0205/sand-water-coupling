#include "StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> StaggeredGridBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec;
	for (int axis = 0; axis < Dim; axis++)
		vec[axis] = _components[axis](pos);
	return vec;
}

template <int Dim>
real StaggeredGridBasedVectorField<Dim>::divergenceAtCellCenter(const VectorDi &cell) const
{
	real acc = 0;
	for (int i = 0; i < StaggeredGrid<Dim>::numberOfCellFaces(); i++) {
		const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
		const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
		const int side = StaggeredGrid<Dim>::cellFaceSide(i);
		acc += _components[axis][face] * side;
	}
	return acc / (2 * _grid->spacing());
}

template <int Dim>
real StaggeredGridBasedVectorField<Dim>::divergence(const VectorDr &pos) const
{
	real div = 0;
	for (const auto [coord, weight] : _grid->cellGrid()->linearIntrplDataPoints(pos))
		div += divergenceAtCellCenter(coord) * weight;
	return div;
}

template class StaggeredGridBasedVectorField<2>;
template class StaggeredGridBasedVectorField<3>;

}
