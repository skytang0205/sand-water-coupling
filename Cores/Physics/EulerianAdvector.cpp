#include "EulerianAdvector.h"

namespace PhysX {

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
	newField.parallelForEach([&](const VectorDi &coord) {
		const VectorDr pos = newField.position(coord);
		newField[coord] = field(backtrace(pos, flow, dt));
	});
	field = newField;
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
	newField.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = newField[axis].position(face);
		newField[axis][face] = field[axis](backtrace(pos, flow, dt));
	});
	field = newField;
}

template <int Dim>
Vector<real, Dim> SemiLagrangianAdvector<Dim>::backtrace(const VectorDr &startPos, const VectorField<Dim> &flow, const real dt) const
{
	const VectorDr midPos = startPos - flow(startPos) * dt * real(0.5);
	VectorDr stopPos = startPos - flow(midPos) * dt;
	return stopPos;
}

template class EulerianAdvector<2>;
template class EulerianAdvector<3>;

template class SemiLagrangianAdvector<2>;
template class SemiLagrangianAdvector<3>;

}
