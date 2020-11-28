#include "GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
real GridBasedScalarField<Dim>::operator()(const VectorDr &pos) const
{
	real val = 0;
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		val += operator[](coord) * weight;
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
	VectorDr grad = VectorDr::Zero();
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		grad += gradientAtDataPoint(coord) * weight;
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
	real lapl = 0;
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		lapl += laplacianAtDataPoint(coord) * weight;
	return lapl;
}

template class GridBasedScalarField<2>;
template class GridBasedScalarField<3>;

}
