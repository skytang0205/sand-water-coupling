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
void EulerianProjector<Dim>::project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights)
{
	buildLinearSystem(velocity, weights);
	solveLinearSystem();
	applyPressureGradient(velocity, weights);
}

template <int Dim>
void EulerianProjector<Dim>::project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &liquidSdf)
{
	buildLinearSystem(velocity, weights, liquidSdf);
	solveLinearSystem();
	applyPressureGradient(velocity, weights, liquidSdf);
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
void EulerianProjector<Dim>::buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights)
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
			const VectorDi face = i & 1 ? cell : nbCell;
			const int nbIdx = int(_reducedPressure.index(nbCell));
			real term = 1;
			if (_reducedPressure.isValid(nbCell)) {
				term *= weights[axis][face];
				diagCoeff += term;
				if (term) _coefficients.push_back(Tripletr(idx, nbIdx, -term));
			}
			div += (i & 1 ? -1 : 1) * term * velocity[axis][face];
		}
		if (!diagCoeff) diagCoeff = 1;
		_velocityDiv[cell] = div;
		_coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (weights[axis][face] > 0) {
			velocity[axis][face] += _reducedPressure[face] - _reducedPressure[face - VectorDi::Unit(axis)];
		}
	});
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &liquidSdf)
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
				const VectorDi face = i & 1 ? cell : nbCell;
				const int nbIdx = int(_reducedPressure.index(nbCell));
				real term = 1;
				if (_reducedPressure.isValid(nbCell)) {
					term *= weights[axis][face];
					if (Surface<Dim>::isInside(liquidSdf[nbCell])) {
						diagCoeff += term;
						if (term) _coefficients.push_back(Tripletr(idx, nbIdx, -term));
					}
					else {
						const real fraction = std::max(Surface<Dim>::theta(liquidSdf[cell], liquidSdf[nbCell]), real(0.01));
						diagCoeff += term / fraction;
					}
				}
				div += (i & 1 ? -1 : 1) * term * velocity[axis][face];
			}
		}
		if (!diagCoeff) diagCoeff = 1;
		_velocityDiv[cell] = div;
		_coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &liquidSdf) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (weights[axis][face] > 0) {
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
