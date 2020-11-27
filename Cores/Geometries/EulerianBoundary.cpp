#include "EulerianBoundary.h"

#include <algorithm>
#include <cstdlib>

namespace PhysX {

template <int Dim>
EulerianBoundary<Dim>::EulerianBoundary(const StaggeredGrid<Dim> *const grid) :
	_domainBox(grid->domainOrigin(), grid->domainLengths()),
	_surface(grid->nodeGrid()),
	_fraction(grid),
	_velocity(grid),
	_normal(grid)
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
		const VectorDr pos = _fraction[axis].position(face);
		if (_fraction.isBoundary(axis, face)) {
			_fraction[axis][face] = 1;
			_velocity[axis][face] = domainBoundaryVelocity ? domainBoundaryVelocity(axis, face) : real(0);
			_normal[axis][face] = -_domainBox.closestNormal(pos)[axis];
		}
		else {
			_fraction[axis][face] = getFaceFraction(axis, face);
			for (const auto &collider : colliders) {
				if (collider->surface()->isInside(pos)) {
					_velocity[axis][face] = collider->velocityAt(pos)[axis];
					_normal[axis][face] = collider->surface()->closestNormal(pos)[axis];
					break;
				}
			}
		}
	});
}

template <int Dim>
void EulerianBoundary<Dim>::enforce(StaggeredGridBasedVectorField<Dim> &fluidVelocity)
{
	fluidVelocity.forEach([&](const int axis, const VectorDi &face) {
		if (_fraction[axis][face] == 1) {
			const VectorDr pos = fluidVelocity[axis].position(face);
			const VectorDr n = _normal(pos).normalized();
			fluidVelocity[axis][face] -= (fluidVelocity(pos) - _velocity(pos)).dot(n) * _normal[axis][face];
		}
	});
}

template <int Dim>
void EulerianBoundary<Dim>::extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, const int maxSteps)
{
	auto newVelocity = fluidVelocity;
	auto visited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(fluidVelocity.staggeredGrid());
	auto newVisited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(fluidVelocity.staggeredGrid());

	visited->parallelForEach([&](const int axis, const VectorDi &face) {
		if (!((*visited)[axis][face] = _fraction[axis][face] < 1))
			newVelocity[axis][face] = 0;
	});

	for (int iter = 0; iter < maxSteps || (maxSteps < 0 && visited->sum<size_t>() < visited->count()); iter++) {
		newVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
			if (!(*visited)[axis][face]) {
				int cnt = 0;
				real sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
					if (newVelocity[axis].isValid(nbFace) && (*visited)[axis][nbFace])
						sum += newVelocity[axis][nbFace], cnt++;
				}
				if (cnt > 0) {
					newVelocity[axis][face] = sum / cnt;
					(*newVisited)[axis][face] = true;
				}
			}
			else (*newVisited)[axis][face] = true;
		});
		visited.swap(newVisited);
	}
	fluidVelocity = newVelocity;
}


template <int Dim>
void EulerianBoundary<Dim>::extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, GridBasedImplicitSurface<Dim> &liquidLevelSet, const int maxSteps)
{
	const auto &liquidSdf = liquidLevelSet.signedDistanceField();
	const auto isLiquidFace = [&](const int axis, const VectorDi &face)->bool {
		const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
		const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
		return _fraction[axis][face] < 1 && (Surface<Dim>::isInside(liquidSdf[cell0]) || Surface<Dim>::isInside(liquidSdf[cell1]));
	};

	auto newVelocity = fluidVelocity;
	newVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (isLiquidFace(axis, face)) return;
		int cnt = 0;
		real sum = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
			if (newVelocity[axis].isValid(nbFace) && isLiquidFace(axis, nbFace))
				sum += newVelocity[axis][nbFace], cnt++;
		}
		if (cnt > 0) newVelocity[axis][face] = sum / cnt;
	});

	const real bandWidth = maxSteps * fluidVelocity.spacing();
	fluidVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = fluidVelocity[axis].position(face);
		if (liquidLevelSet.signedDistance(pos) > 0) {
			fluidVelocity[axis][face] =
				bandWidth < 0 || liquidLevelSet.signedDistance(pos) < bandWidth ? newVelocity[axis](liquidLevelSet.closestPosition(pos)) : real(0);
		}
	});
}

template <int Dim>
real EulerianBoundary<Dim>::getFaceFraction(const int axis, const VectorDi &face) const
{
	const auto &sdf = _surface.signedDistanceField();
	real faceFraction;
	if constexpr (Dim == 2) {
		const real phi0 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 0)];
		const real phi1 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 1)];
		faceFraction = Surface<Dim>::fraction(phi0, phi1);
	}
	else {
		const real phi0 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 0)];
		const real phi1 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 1)];
		const real phi2 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 2)];
		const real phi3 = sdf[StaggeredGrid<Dim>::faceNode(axis, face, 3)];
		const real centerPhi = (phi0 + phi1 + phi2 + phi3) * real(0.25);
		faceFraction = (Surface<Dim>::fraction(centerPhi, phi0, phi1)
			+ Surface<Dim>::fraction(centerPhi, phi1, phi3)
			+ Surface<Dim>::fraction(centerPhi, phi3, phi2)
			+ Surface<Dim>::fraction(centerPhi, phi2, phi0)) * real(0.25);
	}
	return faceFraction > real(0.9) ? real(1) : faceFraction;
}

template class EulerianBoundary<2>;
template class EulerianBoundary<3>;

}
