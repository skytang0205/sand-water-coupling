#pragma once

#include "Geometries/Collider.h"
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

	void project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights);
	void project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi);

protected:

	void solveLinearSystem();

	void buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights);
	void applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights) const;

	void buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi);
	void applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi) const;
};

}
