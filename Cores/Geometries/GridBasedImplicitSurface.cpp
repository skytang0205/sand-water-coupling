#include "GridBasedImplicitSurface.h"

namespace PhysX {

template <int Dim>
void GridBasedImplicitSurface<Dim>::setFromSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &cell) {
		_signedDistanceField[cell] = surface.signedDistance(_signedDistanceField.position(cell));
	});
}

template class GridBasedImplicitSurface<2>;
template class GridBasedImplicitSurface<3>;

}
