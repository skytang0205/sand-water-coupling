#include "GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
real GridBasedScalarField<Dim>::operator()(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += operator[](coords[i]) * weights[i];
	return val;
}

template <int Dim>
Vector<real, Dim> GridBasedScalarField<Dim>::gradientAtDataPoint(const VectorDi &coord) const
{
	VectorDr vec;
	for (int i = 0; i < Dim; i++) {
		vec[i] = operator[](coord[i] < _grid->dataSize()[i] - 1 ? coord + VectorDi::Unit(i) : coord)
			- operator[](coord[i] > 0 ? coord - VectorDi::Unit(i) : coord);
	}
	return vec / (2 * _grid->spacing());
}

template <int Dim>
Vector<real, Dim> GridBasedScalarField<Dim>::gradient(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	VectorDr vec = VectorDr::Zero();
	for (int i = 0; i < (1 << Dim); i++)
		vec += gradientAtDataPoint(coords[i]) * weights[i];
	return vec;
}

template class GridBasedScalarField<2>;
template class GridBasedScalarField<3>;

}
