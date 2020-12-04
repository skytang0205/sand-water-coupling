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
	const LevelSet<Dim> &liquidLevelSet,
	const real surfaceTensionMultiplier)
{
	buildLinearSystem(velocity, boundaryFraction, boundaryVelocity, liquidLevelSet, surfaceTensionMultiplier);
	solveLinearSystem();
	applyPressureGradient(velocity, boundaryFraction, liquidLevelSet, surfaceTensionMultiplier);
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
	_reducedPressure.forEach([&](const VectorDi &cell) {
		const int idx = int(_reducedPressure.index(cell));
		real diagCoeff = 0;
		real div = 0;
		for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
			const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
			const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
			const int side = StaggeredGrid<Dim>::cellFaceSide(i);
			const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
			const real weight = 1 - boundaryFraction[axis][face];
			if (weight > 0) {
				diagCoeff += weight;
				_coefficients.push_back(Tripletr(idx, int(_reducedPressure.index(nbCell)), -weight));
				div += side * weight * velocity[axis][face];
			}
			if (weight < 1)
				div += side * (1 - weight) * boundaryVelocity[axis][face];
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
			const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
			const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
			velocity[axis][face] += _reducedPressure[cell1] - _reducedPressure[cell0];
		}
	});
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(
	StaggeredGridBasedVectorField<Dim> &velocity,
	const StaggeredGridBasedData<Dim> &boundaryFraction,
	const StaggeredGridBasedVectorField<Dim> &boundaryVelocity,
	const LevelSet<Dim> &liquidLevelSet,
	const real surfaceTensionMultiplier)
{
	_coefficients.clear();

	const auto &liquidSdf = liquidLevelSet.signedDistanceField();
	_reducedPressure.forEach([&](const VectorDi &cell) {
		const int idx = int(_reducedPressure.index(cell));
		real diagCoeff = 0;
		real div = 0;
		if (Surface<Dim>::isInside(liquidSdf[cell])) {
			for (int i = 0; i < Grid<Dim>::numberOfNeighbors(); i++) {
				const VectorDi nbCell = Grid<Dim>::neighbor(cell, i);
				const int axis = StaggeredGrid<Dim>::cellFaceAxis(i);
				const int side = StaggeredGrid<Dim>::cellFaceSide(i);
				const VectorDi face = StaggeredGrid<Dim>::cellFace(cell, i);
				const real weight = 1 - boundaryFraction[axis][face];
				if (weight > 0) {
					if (Surface<Dim>::isInside(liquidSdf[nbCell])) {
						diagCoeff += weight;
						_coefficients.push_back(Tripletr(idx, int(_reducedPressure.index(nbCell)), -weight));
					}
					else {
						const real theta = Surface<Dim>::theta(liquidSdf[cell], liquidSdf[nbCell]);
						const real intfCoef = 1 / std::max(theta, real(0.001));
						diagCoeff += weight * intfCoef;
						if (surfaceTensionMultiplier)
							div -= getReducedPressureJump(cell, nbCell, theta, liquidLevelSet, surfaceTensionMultiplier) * weight * intfCoef;
					}
					div += side * weight * velocity[axis][face];
				}
				if (weight < 1)
					div += side * (1 - weight) * boundaryVelocity[axis][face];
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
	const LevelSet<Dim> &liquidLevelSet,
	const real surfaceTensionMultiplier) const
{
	const auto &liquidSdf = liquidLevelSet.signedDistanceField();
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
		if (boundaryFraction[axis][face] < 1) {
			const VectorDi cell0 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 0);
			const VectorDi cell1 = StaggeredGrid<Dim>::faceAdjacentCell(axis, face, 1);
			const real phi0 = liquidSdf[cell0];
			const real phi1 = liquidSdf[cell1];
			if (Surface<Dim>::isInside(phi0) || Surface<Dim>::isInside(phi1)) {
				const real intfCoef = 1 / std::max(Surface<Dim>::fraction(phi0, phi1), real(0.001));
				velocity[axis][face] += (_reducedPressure[cell1] - _reducedPressure[cell0]) * intfCoef;
				if (surfaceTensionMultiplier && Surface<Dim>::isInterface(phi0, phi1)) {
					const real theta = Surface<Dim>::theta(phi0, phi1);
					const real reducedPressureJump = getReducedPressureJump(cell0, cell1, theta, liquidLevelSet, surfaceTensionMultiplier) * intfCoef;
					velocity[axis][face] += (Surface<Dim>::isInside(phi0) ? -1 : 1) * reducedPressureJump;
				}
			}
		}
	});
}

template <int Dim>
real EulerianProjector<Dim>::getReducedPressureJump(
	const VectorDi &cell0,
	const VectorDi &cell1,
	const real theta,
	const LevelSet<Dim> &liquidLevelSet,
	const real surfaceTensionMultiplier) const
{
	const VectorDr pos = (1 - theta) * _reducedPressure.position(cell0) + theta * _reducedPressure.position(cell1);
	return liquidLevelSet.curvature(pos) * surfaceTensionMultiplier;
}

template class EulerianProjector<2>;
template class EulerianProjector<3>;

}
