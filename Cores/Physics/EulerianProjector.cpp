#include "EulerianProjector.h"

namespace PhysX {

template <int Dim>
EulerianProjector<Dim>::EulerianProjector(const Grid<Dim> *const grid) :
	_reducedPressure(grid),
	_velocityDiv(grid),
	_matLaplacian(grid->dataCount(), grid->dataCount()),
	_solver(std::make_unique<IcPCgSolver>())
{ }

template <int Dim>
void EulerianProjector<Dim>::project(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const StaggeredGridBasedVectorField<Dim> &boundaryVelocity)
{
	buildLinearSystem(velocity, boundaryFraction, boundaryVelocity);
	solveLinearSystem();
	applyPressureGradient(velocity, boundaryFraction);
}

template <int Dim>
void EulerianProjector<Dim>::project(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const StaggeredGridBasedVectorField<Dim> &boundaryVelocity,
	const GridBasedScalarField<Dim> &liquidSdf)
{
	buildLinearSystem(velocity, boundaryFraction, boundaryVelocity, liquidSdf);
	solveLinearSystem();
	applyPressureGradient(velocity, boundaryFraction, liquidSdf);
}

template <int Dim>
void EulerianProjector<Dim>::solveLinearSystem()
{
	_solver->solve(
		_matLaplacian,
		Eigen::Map<VectorXr, Eigen::Aligned>(_reducedPressure.data(), _reducedPressure.count()),
		Eigen::Map<VectorXr, Eigen::Aligned>(_velocityDiv.data(), _velocityDiv.count()));
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const StaggeredGridBasedVectorField<Dim> &boundaryVelocity)
{
	_coefficients.clear();
	_matLaplacian.setZero();
	_reducedPressure.forEach([&](const VectorDi &cell) {
		const int idx = int(_reducedPressure.index(cell));
		real diagCoeff = 0;
		real div = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = i >> 1;
			const int dir = i & 1 ? -1 : 1;
			const VectorDi face = dir < 0 ? cell : nbCell;
			const real weight = 1 - boundaryFraction[axis][face];
			if (weight > 0) {
				diagCoeff += weight;
				_coefficients.push_back(Tripletr(idx, int(_reducedPressure.index(nbCell)), -weight));
				div += dir * weight * velocity[axis][face];
			}
			if (weight < 1)
				div += dir * (1 - weight) * boundaryVelocity[axis][face];
		}
		if (!diagCoeff) diagCoeff = 1;
		_velocityDiv[cell] = div;
		_coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (boundaryFraction[axis][face] < 1) {
			velocity[axis][face] += _reducedPressure[face] - _reducedPressure[face - VectorDi::Unit(axis)];
		}
	});
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const StaggeredGridBasedVectorField<Dim> &boundaryVelocity,
	const GridBasedScalarField<Dim> &liquidSdf)
{
	_coefficients.clear();
	_matLaplacian.setZero();
	_reducedPressure.forEach([&](const VectorDi &cell) {
		const int idx = int(_reducedPressure.index(cell));
		real diagCoeff = 0;
		real div = 0;
		if (Surface<Dim>::isInside(liquidSdf[cell])) {
			for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
				const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
				const int axis = i >> 1;
				const int dir = i & 1 ? -1 : 1;
				const VectorDi face = dir < 0 ? cell : nbCell;
				const real weight = 1 - boundaryFraction[axis][face];
				if (weight > 0) {
					if (Surface<Dim>::isInside(liquidSdf[nbCell])) {
						diagCoeff += weight;
						_coefficients.push_back(Tripletr(idx, int(_reducedPressure.index(nbCell)), -weight));
					}
					else {
						const real theta = std::max(Surface<Dim>::theta(liquidSdf[cell], liquidSdf[nbCell]), real(0.01));
						diagCoeff += weight / theta;
					}
					div += dir * weight * velocity[axis][face];
				}
				if (weight < 1)
					div += dir * (1 - weight) * boundaryVelocity[axis][face];
			}
		}
		if (!diagCoeff) diagCoeff = 1;
		_velocityDiv[cell] = div;
		_coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const GridBasedScalarField<Dim> &liquidSdf) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (boundaryFraction[axis][face] < 1) {
			const VectorDi cell0 = face - VectorDi::Unit(axis);
			const VectorDi cell1 = face;
			if (Surface<Dim>::isInside(liquidSdf[cell0]) || Surface<Dim>::isInside(liquidSdf[cell1])) {
				const real fraction = std::max(Surface<Dim>::fraction(liquidSdf[cell0], liquidSdf[cell1]), real(0.01));
				velocity[axis][face] += (_reducedPressure[cell1] - _reducedPressure[cell0]) / fraction;
			}
		}
	});
}

template class EulerianProjector<2>;
template class EulerianProjector<3>;

}
