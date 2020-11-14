#pragma once

#include "Geometries/Collider.h"
#include "Structures/GridBasedScalarField.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class EulerianProjector
{
	DECLARE_DIM_TYPES(Dim)

protected:

	GridBasedScalarField<Dim> _reducedPressure; // reducedPressure = dt / dx / rho * pressure

public:

	EulerianProjector(const Grid<Dim> *const grid) : _reducedPressure(grid) { }

	EulerianProjector(const EulerianProjector &rhs) = delete;
	EulerianProjector &operator=(const EulerianProjector &rhs) = delete;
	virtual ~EulerianProjector() = default;

	void project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights);

protected:

	void buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights);
	void applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights) const;
};

}
