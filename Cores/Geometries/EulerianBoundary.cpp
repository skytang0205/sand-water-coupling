#include "EulerianBoundary.h"

namespace PhysX {

template <int Dim>
EulerianBoundary<Dim>::EulerianBoundary(const StaggeredGrid<Dim> *const grid) :
	_surface(grid->nodeGrid()),
	_fraction(grid),
	_velocity(grid)
{ }

template <int Dim>
void EulerianBoundary<Dim>::reset(
	const std::vector<std::unique_ptr<Collider<Dim>>> &colliders,
	const std::function<real(const int axis, const VectorDi &face)> &domainBoundaryVelocity)
{
	for (const auto &collider : colliders)
		_surface.unionSurface(*collider->surface());
	// Colliders must not intersect each other, or we should reinitialize _surface here.

	_fraction.parallelForEach([&](const int axis, const VectorDi &face) {
		if (_fraction.isBoundary(axis, face)) {
			_fraction[axis][face] = 1;
			_velocity[axis][face] = domainBoundaryVelocity ? domainBoundaryVelocity(axis, face) : real(0);
		}
		else {
			_fraction[axis][face] = getFaceFraction(axis, face);
			_velocity[axis][face] = 0;
			const VectorDr pos = _fraction[axis].position(face);
			for (const auto &collider : colliders) {
				if (collider->surface()->isInside(pos)) {
					_velocity[axis][face] = collider->velocityAt(pos)[axis];
					break;
				}
			}
		}
	});
}

template class EulerianBoundary<2>;
template class EulerianBoundary<3>;

}
