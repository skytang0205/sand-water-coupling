#include "EulerianProjector.h"

namespace PhysX {

template <int Dim>
void EulerianProjector<Dim>::project(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights)
{
	buildLinearSystem(velocity, weights);
	// TODO: solve the linear system
	applyPressureGradient(velocity, weights);
}

template <int Dim>
void EulerianProjector<Dim>::buildLinearSystem(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights)
{
}

template <int Dim>
void EulerianProjector<Dim>::applyPressureGradient(StaggeredGridBasedVectorField<Dim> &velocity, const StaggeredGridBasedData<Dim> &weights) const
{
	velocity.parallelForEach([&](const int axis, const VectorDi &face) {
			if (!velocity.isBoundary(axis, face) && weights[axis][face] > 0) {
				velocity[axis][face] -= _reducedPressure[face] - _reducedPressure[face - VectorDi::Unit(axis)];
			}
		});
}

}
