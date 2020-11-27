#include "EulerianAdvector.h"

namespace PhysX {

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	GridBasedScalarField<Dim> newField(field.grid());
	advect(field, newField, flow, dt);
	field = newField;
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	StaggeredGridBasedVectorField<Dim> newField(field.staggeredGrid());
	advect(field, newField, flow, dt);
	field = newField;
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(ParticlesVectorAttribute<Dim> &positions, const VectorField<Dim> &flow, const real dt)
{
	positions.parallelForEach([&](const int i) {
		positions[i] = trace(positions[i], flow, dt);
	});
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(const GridBasedScalarField<Dim> &field, GridBasedScalarField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const
{
	newField.parallelForEach([&](const VectorDi &coord) {
		const VectorDr pos = newField.position(coord);
		newField[coord] = field(trace(pos, flow, -dt));
	});
}

template <int Dim>
void SemiLagrangianAdvector<Dim>::advect(const StaggeredGridBasedVectorField<Dim> &field, StaggeredGridBasedVectorField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const
{
	newField.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = newField[axis].position(face);
		newField[axis][face] = field[axis](trace(pos, flow, -dt));
	});
}

template <int Dim>
Vector<Dim, real> SemiLagrangianAdvector<Dim>::trace(const VectorDr &startPos, const VectorField<Dim> &flow, const real dt) const
{
	const VectorDr midPos = startPos + flow(startPos) * dt * real(0.5);
	VectorDr stopPos = startPos + flow(midPos) * dt;
	return stopPos;
}

template <int Dim>
void MacCormackAdvector<Dim>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	GridBasedScalarField<Dim> forwardField(field.grid());
	GridBasedScalarField<Dim> backwardField(field.grid());
	SemiLagrangianAdvector<Dim>::advect(field, forwardField, flow, dt);
	SemiLagrangianAdvector<Dim>::advect(forwardField, backwardField, flow, -dt);
	field.parallelForEach([&](const VectorDi &coord) {
		field[coord] = forwardField[coord] + (field[coord] - backwardField[coord]) * real(0.5);
	});
}

template <int Dim>
void MacCormackAdvector<Dim>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	StaggeredGridBasedVectorField<Dim> forwardField(field.staggeredGrid());
	StaggeredGridBasedVectorField<Dim> backwardField(field.staggeredGrid());
	SemiLagrangianAdvector<Dim>::advect(field, forwardField, flow, dt);
	SemiLagrangianAdvector<Dim>::advect(forwardField, backwardField, flow, -dt);
	field.parallelForEach([&](const int axis, const VectorDi &coord) {
		field[axis][coord] = forwardField[axis][coord] + (field[axis][coord] - backwardField[axis][coord]) * real(0.5);
	});
}

template <int Dim>
void MacCormackAdvector<Dim>::advect(ParticlesVectorAttribute<Dim> &positions, const VectorField<Dim> &flow, const real dt)
{
	positions.parallelForEach([&](const int i) {
		const VectorDr forwardPos = trace(positions[i], flow, dt);
		const VectorDr backwardPos = trace(forwardPos, flow, -dt);
		positions[i] = forwardPos + (positions[i] - backwardPos) * real(0.5);
	});
}

template class EulerianAdvector<2>;
template class EulerianAdvector<3>;

template class SemiLagrangianAdvector<2>;
template class SemiLagrangianAdvector<3>;

template class MacCormackAdvector<2>;
template class MacCormackAdvector<3>;

}
