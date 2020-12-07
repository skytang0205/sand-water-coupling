#include "GridBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> GridBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec = VectorDr::Zero();
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		vec += at(coord) * weight;
	return vec;
}

template <int Dim>
real GridBasedVectorField<Dim>::divergenceAtDataPoint(const VectorDi &coord) const
{
	real acc = 0;
	for (int i = 0; i < Dim; i++) {
		acc += at(coord + VectorDi::Unit(i))[i] - at(coord - VectorDi::Unit(i))[i];
	}
	return acc * real(.5) * _grid->invSpacing();
}

template <int Dim>
real GridBasedVectorField<Dim>::divergence(const VectorDr &pos) const
{
	real div = 0;
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		div += divergenceAtDataPoint(coord) * weight;
	return div;
}

template class GridBasedVectorField<2>;
template class GridBasedVectorField<3>;


}
