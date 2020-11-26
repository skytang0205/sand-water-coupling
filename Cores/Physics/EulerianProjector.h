#pragma once

#include "Geometries/Collider.h"
#include "Geometries/GridBasedImplicitSurface.h"
#include "Solvers/SparseSolver.h"
#include "Structures/GridBasedScalarField.h"
#include "Structures/StaggeredGridBasedData.h"
#include "Structures/StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class EulerianProjector
{
	DECLARE_DIM_TYPES(Dim)

protected:

	GridBasedScalarField<Dim> _reducedPressure; // reducedPressure = -dt / dx / rho * pressure
	GridBasedScalarField<Dim> _velocityDiv; // velocityDiv = Div(velocity) * dx

	std::vector<Tripletr> _coefficients;
	SparseMatrixr _matLaplacian;

	std::unique_ptr<SparseSolver> _solver;

public:

	EulerianProjector(const Grid<Dim> *const grid);

	EulerianProjector(const EulerianProjector &rhs) = delete;
	EulerianProjector &operator=(const EulerianProjector &rhs) = delete;
	virtual ~EulerianProjector() = default;

	void project(
		StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction,
		const StaggeredGridBasedVectorField<Dim> &boundaryVelocity);
	void project(StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction,
		const StaggeredGridBasedVectorField<Dim> &boundaryVelocity,
		const GridBasedImplicitSurface<Dim> &liquidSurface,
		const real surfaceTensionMultiplier = 0);

protected:

	void solveLinearSystem();

	void buildLinearSystem(
		StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction,
		const StaggeredGridBasedVectorField<Dim> &boundaryVelocity);
	void applyPressureGradient(
		StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction) const;

	void buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction,
		const StaggeredGridBasedVectorField<Dim> &boundaryVelocity,
		const GridBasedImplicitSurface<Dim> &liquidSurface,
		const real surfaceTensionMultiplier);
	void applyPressureGradient(
		StaggeredGridBasedVectorField<Dim> &velocity,
		const StaggeredGridBasedData<Dim> &boundaryFraction,
		const GridBasedImplicitSurface<Dim> &liquidSurface,
		const real surfaceTensionMultiplier) const;

	real getReducedPressureJump(
		const VectorDi &cell0,
		const VectorDi &cell1,
		const real theta,
		const GridBasedImplicitSurface<Dim> &liquidSurface,
		const real surfaceTensionMultiplier) const;
};

}
