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
Vector<Dim, real> GridBasedScalarField<Dim>::gradientAtDataPoint(const VectorDi &coord) const
{
	VectorDr acc;
	for (int i = 0; i < Dim; i++) {
		acc[i] = operator[](coord[i] < _grid->dataSize()[i] - 1 ? coord + VectorDi::Unit(i) : coord)
			- operator[](coord[i] > 0 ? coord - VectorDi::Unit(i) : coord);
	}
	return acc / (2 * _grid->spacing());
}

template <int Dim>
Vector<Dim, real> GridBasedScalarField<Dim>::gradient(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	VectorDr grad = VectorDr::Zero();
	for (int i = 0; i < (1 << Dim); i++)
		grad += gradientAtDataPoint(coords[i]) * weights[i];
	return grad;
}

template <int Dim>
real GridBasedScalarField<Dim>::laplacianAtDataPoint(const VectorDi &coord) const
{
	real acc = 0;
	const real centerVal = operator[](coord);
	for (int i = 0; i < Dim; i++) {
		if (coord[i] < _grid->dataSize()[i] - 1)
			acc += operator[](coord + VectorDi::Unit(i)) - centerVal;
		if (coord[i] > 0)
			acc -= centerVal - operator[](coord - VectorDi::Unit(i));
	}
	return acc / (_grid->spacing() * _grid->spacing());
}

template <int Dim>
real GridBasedScalarField<Dim>::laplacian(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	real lapl = 0;
	for (int i = 0; i < (1 << Dim); i++)
		lapl += laplacianAtDataPoint(coords[i]) * weights[i];
	return lapl;
}

template class GridBasedScalarField<2>;
template class GridBasedScalarField<3>;

}
