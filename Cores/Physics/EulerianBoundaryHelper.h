#pragma once

#include "Geometries/Collider.h"
#include "Geometries/LevelSet.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"
#include "Structures/ParticlesAttribute.h"

#include <functional>
#include <memory>

namespace PhysX {

template <int Dim>
class EulerianBoundaryHelper
{
	DECLARE_DIM_TYPES(Dim)

protected:

	ImplicitBox<Dim> _domainBox;
	LevelSet<Dim> _surface;
	StaggeredGridBasedData<Dim> _fraction;
	StaggeredGridBasedVectorField<Dim> _velocity;
	StaggeredGridBasedVectorField<Dim> _normal;

public:

	EulerianBoundaryHelper(const StaggeredGrid<Dim> *const grid);

	EulerianBoundaryHelper(const EulerianBoundaryHelper &rhs) = delete;
	EulerianBoundaryHelper &operator=(const EulerianBoundaryHelper &rhs) = delete;
	virtual ~EulerianBoundaryHelper() = default;

	const StaggeredGridBasedData<Dim> &fraction() const { return _fraction; }
	const StaggeredGridBasedVectorField<Dim> &velocity() const { return _velocity; }

	void reset(
		const std::vector<std::unique_ptr<Collider<Dim>>> &colliders,
		const std::function<real(const int axis, const VectorDi &face)> &domainBoundaryVelocity);

	void enforce(StaggeredGridBasedVectorField<Dim> &fluidVelocity) const;
	void enforce(ParticlesVectorAttribute<Dim> &markerPositions) const;

	void extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, const int maxSteps = -1) const;
	void extrapolate(StaggeredGridBasedVectorField<Dim> &fluidVelocity, LevelSet<Dim> &liquidLevelSet, const int maxSteps = -1) const;

protected:

	real getFaceFraction(const int axis, const VectorDi &face) const;
};

}
