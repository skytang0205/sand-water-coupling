#include "GridBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> GridBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	VectorDr vec = VectorDr::Zero();
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		vec += operator[](coord) * weight;
	return vec;
}

template <int Dim>
real GridBasedVectorField<Dim>::divergenceAtDataPoint(const VectorDi &coord) const
{
	real acc = 0;
	for (int i = 0; i < Dim; i++) {
		acc += operator[](coord[i] < _grid->dataSize()[i] - 1 ? coord + VectorDi::Unit(i) : coord)[i]
			- operator[](coord[i] > 0 ? coord - VectorDi::Unit(i) : coord)[i];
	}
	return acc / (2 * _grid->spacing());
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
