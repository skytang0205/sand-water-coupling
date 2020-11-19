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
void EulerianProjector<Dim>::project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi)
{
	buildLinearSystem(velocity, weights, phi);
	solveLinearSystem();
	applyPressureGradient(velocity, weights, phi);
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
		for (int axis = 0; axis < Dim; axis++) {
			for (int i = -1; i <= 1; i += 2) {
				const VectorDi nbCell = cell + VectorDi::Unit(axis) * i;
				const VectorDi face = i < 0 ? cell : nbCell;
				const int nbIdx = int(_reducedPressure.index(nbCell));
				real term = 1;
				if (_reducedPressure.isValid(nbCell)) {
					term *= weights[axis][face];
					diagCoeff += term;
					if (term) _coefficients.push_back(Tripletr(idx, nbIdx, -term));
				}
				div += i * term * velocity[axis][face];
			}
		}
		_velocityDiv[cell] = div;
		if (diagCoeff) _coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (!velocity.isBoundary(axis, face) && weights[axis][face] > 0) {
			velocity[axis][face] += _reducedPressure[face] - _reducedPressure[face - VectorDi::Unit(axis)];
		}
	});
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi)
{
	_coefficients.clear();
	_matLaplacian.setZero();
	_reducedPressure.forEach([&](const VectorDi &cell) {
		const int idx = int(_reducedPressure.index(cell));
		real diagCoeff = 0;
		real div = 0;
		if (Surface<Dim>::inside(phi[cell])) {
			for (int axis = 0; axis < Dim; axis++) {
				for (int i = -1; i <= 1; i += 2) {
					const VectorDi nbCell = cell + VectorDi::Unit(axis) * i;
					const VectorDi face = i < 0 ? cell : nbCell;
					const int nbIdx = int(_reducedPressure.index(nbCell));
					real term = 1;
					if (_reducedPressure.isValid(nbCell)) {
						term *= weights[axis][face];
						if (Surface<Dim>::inside(phi[nbCell])) {
							diagCoeff += term;
							if (term) _coefficients.push_back(Tripletr(idx, nbIdx, -term));
						}
						else {
							const real fraction = std::max(Surface<Dim>::theta(phi[cell], phi[nbCell]), real(0.01));
							diagCoeff += term / fraction;
						}
					}
					div += i * term * velocity[axis][face];
				}
			}
		}
		else diagCoeff = 1;
		_velocityDiv[cell] = div;
		if (diagCoeff) _coefficients.push_back(Tripletr(idx, idx, diagCoeff));
	});
	_matLaplacian.setFromTriplets(_coefficients.begin(), _coefficients.end());
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights, const GridBasedScalarField<Dim> &phi) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (!velocity.isBoundary(axis, face) && weights[axis][face] > 0) {
			const VectorDi cell0 = face - VectorDi::Unit(axis);
			const VectorDi cell1 = face;
			if (Surface<Dim>::inside(phi[cell0]) || Surface<Dim>::inside(phi[cell1])) {
				const real fraction = std::max(Surface<Dim>::fraction(phi[cell0], phi[cell1]), real(0.01));
				velocity[axis][face] += (_reducedPressure[face] - _reducedPressure[face - VectorDi::Unit(axis)]) / fraction;
			}
		}
	});
}

template class EulerianProjector<2>;
template class EulerianProjector<3>;

}
