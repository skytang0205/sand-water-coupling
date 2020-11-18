#include "EulerianAdvector.h"

namespace PhysX {

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
	newField.parallelForEach([&](const VectorDi &coord) {
		const VectorDr startPos = newField.position(coord);
		const VectorDr midPos = startPos - flow(startPos) * dt * real(0.5);
		newField[coord] = field(startPos - flow(midPos) * dt);
	});
	field = newField;
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
	newField.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr startPos = newField[axis].position(face);
		const VectorDr midPos = startPos - flow(startPos) * dt * real(0.5);
		newField[axis][face] = field[axis](startPos - flow(midPos) * dt);
	});
	field = newField;
}

template class EulerianAdvector<2>;
template class EulerianAdvector<3>;

template class SemiLagrangianAdvector<2>;
template class SemiLagrangianAdvector<3>;

}
