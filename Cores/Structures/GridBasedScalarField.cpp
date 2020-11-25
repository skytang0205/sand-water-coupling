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
	VectorDr val;
	for (int i = 0; i < Dim; i++) {
		val[i] = operator[](coord[i] < _grid->dataSize()[i] - 1 ? coord + VectorDi::Unit(i) : coord)
			- operator[](coord[i] > 0 ? coord - VectorDi::Unit(i) : coord);
	}
	return val / (2 * _grid->spacing());
}

template <int Dim>
Vector<real, Dim> GridBasedScalarField<Dim>::gradient(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	VectorDr val = VectorDr::Zero();
	for (int i = 0; i < (1 << Dim); i++)
		val += gradientAtDataPoint(coords[i]) * weights[i];
	return val;
}

template <int Dim>
real GridBasedScalarField<Dim>::laplacianAtDataPoint(const VectorDi &coord) const
{
	real val = 0;
	for (int i = 0; i < Dim; i++) {
		val += (operator[](coord[i] < _grid->dataSize()[i] - 1 ? coord + VectorDi::Unit(i) : coord) - operator[](coord))
			- (operator[](coord) - operator[](coord[i] > 0 ? coord - VectorDi::Unit(i) : coord));
	}
	return val / (_grid->spacing() * _grid->spacing());
}

template <int Dim>
real GridBasedScalarField<Dim>::laplacian(const VectorDr &pos) const
{
	std::array<VectorDi, 1 << Dim> coords;
	std::array<real, 1 << Dim> weights;
	_grid->getLerpCoordsAndWeights(pos, coords, weights);

	real val = 0;
	for (int i = 0; i < (1 << Dim); i++)
		val += laplacianAtDataPoint(coords[i]) * weights[i];
	return val;
}

template class GridBasedScalarField<2>;
template class GridBasedScalarField<3>;

}
