#include "GridBasedVectorField.h"

namespace PhysX {

template <int Dim>
Vector<Dim, real> GridBasedVectorField<Dim>::operator()(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	VectorDr vec = VectorDr::Zero();
	for (int i = 0; i < (1 << Dim); i++)
		vec += operator[](coords[i]) * weights[i];
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
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	real div = 0;
	for (int i = 0; i < (1 << Dim); i++)
		div += divergenceAtDataPoint(coords[i]) * weights[i];
	return div;
}

template class GridBasedVectorField<2>;
template class GridBasedVectorField<3>;


}
