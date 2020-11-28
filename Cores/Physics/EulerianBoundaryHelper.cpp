#include "EulerianBoundaryHelper.h"

#include <algorithm>
#include <cstdlib>

namespace PhysX {

template <int Dim>
EulerianBoundaryHelper<Dim>::EulerianBoundaryHelper(const StaggeredGrid<Dim> *const grid) :
	_domainBox(grid->domainOrigin(), grid->domainLengths()),
	_surface(grid->nodeGrid()),
	_fraction(grid),
	_velocity(grid),
	_normal(grid)
{ }

template <int Dim>
void EulerianBoundaryHelper<Dim>::reset(
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
void EulerianBoundaryHelper<Dim>::enforce(StaggeredGridBasedVectorField<Dim> &fluidVelocity) const
{
	auto newFluidVelocity = fluidVelocity;
	newFluidVelocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (_fraction[axis][face] == 1) {
			const VectorDr pos = newFluidVelocity[axis].position(face);
			const VectorDr n = _normal(pos).normalized();
			if (n.any())
				newFluidVelocity[axis][face] -= (fluidVelocity(pos) - _velocity(pos)).dot(n) * _normal[axis][face];
			else
				newFluidVelocity[axis][face] = _velocity[axis][face];
		}
	});
	fluidVelocity = newFluidVelocity;
}

template <int Dim>
void EulerianBoundaryHelper<Dim>::enforce(ParticlesVectorAttribute<Dim> &markerPositions) const
{
	markerPositions.parallelForEach([&](const int i) {
		VectorDr &pos = markerPositions[i];
		if (!_domainBox.isInside(pos))
			pos = _domainBox.closestPosition(pos);
		else if (_surface.isInside(pos))
			pos = _surface.closestPosition(pos);
	});
}

template <int Dim>
void EulerianBoundaryHelper<Dim>::extrapolate(StaggeredGridBasedVectorField<Dim> &velocity, const int maxSteps) const
{
	auto newVelocty = velocity;
	auto visited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(velocity.staggeredGrid());
	auto newVisited = std::make_unique<StaggeredGridBasedData<Dim, uchar>>(velocity.staggeredGrid());

	visited->parallelForEach([&](const int axis, const VectorDi &face) {
		if (!((*visited)[axis][face] = _fraction[axis][face] < 1))
			newVelocty[axis][face] = 0;
	});

	for (int iter = 0; iter < maxSteps || (maxSteps < 0 && visited->sum<size_t>() < visited->count()); iter++) {
		newVelocty.parallelForEach([&](const int axis, const VectorDi &face) {
			if (!(*visited)[axis][face]) {
				int cnt = 0;
				real sum = 0;
				for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
					const VectorDi &nbFace = Grid<Dim>::neighbor(face, i);
					if (newVelocty[axis].isValid(nbFace) && (*visited)[axis][nbFace])
						sum += newVelocty[axis][nbFace], cnt++;
				}
				if (cnt > 0) {
					newVelocty[axis][face] = sum / cnt;
					(*newVisited)[axis][face] = true;
				}
			}
			else (*newVisited)[axis][face] = true;
		});
		visited.swap(newVisited);
	}
	velocity = newVelocty;
}


template <int Dim>
void EulerianBoundaryHelper<Dim>::extrapolate(StaggeredGridBasedVectorField<Dim> &velocity, LevelSet<Dim> &liquidLevelSet, const int maxSteps) const
{
	const auto &liquidSdf = liquidLevelSet.signedDistanceField();
	const auto isLiquidFace = [&](const int axis, const VectorDi &face)->bool {
		const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
		const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
		return _fraction[axis][face] < 1 && (Surface<Dim>::isInside(liquidSdf[cell0]) || Surface<Dim>::isInside(liquidSdf[cell1]));
	};

	auto newVelocity = velocity;
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

	const real bandWidth = maxSteps * velocity.spacing();
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = velocity[axis].position(face);
		if (!isLiquidFace(axis, face)) {
			velocity[axis][face] =
				bandWidth < 0 || liquidLevelSet.signedDistance(pos) < bandWidth ? newVelocity[axis](liquidLevelSet.closestPosition(pos)) : real(0);
		}
	});
}

template <int Dim>
real EulerianBoundaryHelper<Dim>::getFaceFraction(const int axis, const VectorDi &face) const
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

template class EulerianBoundaryHelper<2>;
template class EulerianBoundaryHelper<3>;

}
