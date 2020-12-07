#include "GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
real GridBasedScalarField<Dim>::operator()(const VectorDr &pos) const
{
	real val = 0;
	for (const auto [coord, weight] : _grid->linearIntrplDataPoints(pos))
		val += at(coord) * weight;
	return val;
}

template <int Dim>
Vector<Dim, real> GridBasedScalarField<Dim>::gradientAtDataPoint(const VectorDi &coord) const
{
	VectorDr acc;
	for (int i = 0; i < Dim; i++) {
		acc[i] = at(coord + VectorDi::Unit(i)) - at(coord - VectorDi::Unit(i));
	}
	return acc * real(.5) * _grid->invSpacing();
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
	const real centerVal = at(coord);
	for (int i = 0; i < Dim; i++) {
		acc += at(coord + VectorDi::Unit(i)) - centerVal;
		acc -= centerVal - at(coord - VectorDi::Unit(i));
	}
	return acc * _grid->invSpacing() * _grid->invSpacing();
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
