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

template class GridBasedImplicitSurface<2>;
template class GridBasedImplicitSurface<3>;

}
