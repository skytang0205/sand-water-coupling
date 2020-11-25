#pragma once

#include "Geometries/Collider.h"
#include "Geometries/GridBasedImplicitSurface.h"
#include "Geometries/LevelSetReinitializer.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"

#include <functional>
#include <memory>

namespace PhysX {

template <int Dim>
class EulerianBoundary
{
	DECLARE_DIM_TYPES(Dim)

protected:

	ImplicitBox<Dim> _domainBox;
	GridBasedImplicitSurface<Dim> _surface;
	StaggeredGridBasedData<Dim> _fraction;
	StaggeredGridBasedVectorField<Dim> _velocity;
	StaggeredGridBasedVectorField<Dim> _normal;

public:

	EulerianBoundary(const StaggeredGrid<Dim> *const grid);

	EulerianBoundary(const EulerianBoundary &rhs) = delete;
	EulerianBoundary &operator=(const EulerianBoundary &rhs) = delete;
	virtual ~EulerianBoundary() = default;

	const StaggeredGridBasedData<Dim> &fraction() const { return _fraction; }
	const StaggeredGridBasedVectorField<Dim> &velocity() const { return _velocity; }

	void reset(
		const std::vector<std::unique_ptr<Collider<Dim>>> &colliders,
		const std::function<real(const int axis, const VectorDi &face)> &domainBoundaryVelocity);

	void enforce(StaggeredGridBasedVectorField<Dim> &fluidVelocity);

	void extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, const int maxIterations = -1);
	void extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, GridBasedImplicitSurface<Dim> &liquidLevelSet, const int maxIterations = -1);

protected:

	real getFaceFraction(const int axis, const VectorDi &face) const;
};

}
