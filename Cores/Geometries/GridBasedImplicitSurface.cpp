#include "GridBasedImplicitSurface.h"

namespace PhysX {

template <int Dim>
void GridBasedImplicitSurface<Dim>::setFromSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = surface.signedDistance(_signedDistanceField.position(coord));
	});
}

template class GridBasedImplicitSurface<2>;
template class GridBasedImplicitSurface<3>;

}
