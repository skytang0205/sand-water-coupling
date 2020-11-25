#include "GridBasedImplicitSurface.h"

namespace PhysX {

template <int Dim>
void GridBasedImplicitSurface<Dim>::unionSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::min(_signedDistanceField[coord], surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
void GridBasedImplicitSurface<Dim>::intersectSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::max(_signedDistanceField[coord], surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
void GridBasedImplicitSurface<Dim>::exceptSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::max(_signedDistanceField[coord], -surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
real GridBasedImplicitSurface<Dim>::curvature(const VectorDr &pos) const
{
	const real dx = _signedDistanceField.spacing();
	for (int i = 0; i < Dim; i++) {
		acc += closestNormal(pos + VectorDr::Unit(i) * dx / 2)[i] - closestNormal(pos - VectorDr::Unit(i) * dx / 2)[i];
	}
	acc /= dx;
	return std::abs(acc) < 1 / dx ? acc : (acc < 0 ? -1 : 1) / dx;
}

template class GridBasedImplicitSurface<2>;
template class GridBasedImplicitSurface<3>;

}
