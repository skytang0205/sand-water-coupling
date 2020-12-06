#include "EulerianAdvector.h"

namespace PhysX {

template <int Dim, int RungeKuttaOrder>
void SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	GridBasedScalarField<Dim> newField(field.grid());
	advect(field, newField, flow, dt);
	field = newField;
}

template <int Dim, int RungeKuttaOrder>
void SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	StaggeredGridBasedVectorField<Dim> newField(field.staggeredGrid());
	advect(field, newField, flow, dt);
	field = newField;
}

template <int Dim, int RungeKuttaOrder>
void SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(Particles<Dim> &particles, const VectorField<Dim> &flow, const real dt)
{
	particles.parallelForEach([&](const int i) {
		VectorDr &pos = particles.positions[i];
		pos = trace(pos, flow, dt);
	});
}

template <int Dim, int RungeKuttaOrder>
void SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(const GridBasedScalarField<Dim> &field, GridBasedScalarField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const
{
	newField.parallelForEach([&](const VectorDi &coord) {
		const VectorDr pos = newField.position(coord);
		newField[coord] = field(trace(pos, flow, -dt));
	});
}

template <int Dim, int RungeKuttaOrder>
void SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(const StaggeredGridBasedVectorField<Dim> &field, StaggeredGridBasedVectorField<Dim> &newField, const VectorField<Dim> &flow, const real dt) const
{
	newField.parallelForEach([&](const int axis, const VectorDi &face) {
		const VectorDr pos = newField[axis].position(face);
		newField[axis][face] = field[axis](trace(pos, flow, -dt));
	});
}

template <int Dim, int RungeKuttaOrder>
Vector<Dim, real> SemiLagrangianAdvector<Dim, RungeKuttaOrder>::trace(const VectorDr &startPos, const VectorField<Dim> &flow, const real dt) const
{
	if constexpr (RungeKuttaOrder == 1) return startPos + flow(startPos) * dt;
	else if constexpr (RungeKuttaOrder == 2) {
		const VectorDr vel0 = flow(startPos);
		const VectorDr vel1 = flow(startPos + vel0 * dt * 2 / 3);
		return startPos + (vel0 + 3 * vel1) * dt / 4;
	}
	else if constexpr (RungeKuttaOrder == 3) {
		const VectorDr vel0 = flow(startPos);
		const VectorDr vel1 = flow(startPos + vel0 * dt / 2);
		const VectorDr vel2 = flow(startPos + vel1 * dt * 3 / 4);
		return startPos + (vel0 * 2 + vel1 * 3 + vel2 * 4) * dt / 9;
	}
	else {
		const VectorDr vel0 = flow(startPos);
		const VectorDr vel1 = flow(startPos + vel0 * dt * real(.4));
		const VectorDr vel2 = flow(startPos + (vel0 * real(.29697761) + vel1 * real(.15875964)) * dt);
		const VectorDr vel3 = flow(startPos + (vel0 * real(.21810040) - vel1 * real(3.05096516) + vel2 * (3.83286476)) * dt);
		return startPos + (vel0 * real(.17476028) - vel1 * (.55148066) + vel2 * real(1.20553560) + vel3 * real(.17118478)) * dt;
	}
}

template <int Dim, int RungeKuttaOrder>
void MacCormackAdvector<Dim, RungeKuttaOrder>::advect(GridBasedScalarField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	GridBasedScalarField<Dim> forwardField(field.grid());
	GridBasedScalarField<Dim> backwardField(field.grid());
	SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(field, forwardField, flow, dt);
	SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(forwardField, backwardField, flow, -dt);
	field.parallelForEach([&](const VectorDi &coord) {
		field[coord] = forwardField[coord] + (field[coord] - backwardField[coord]) * real(0.5);
	});
}

template <int Dim, int RungeKuttaOrder>
void MacCormackAdvector<Dim, RungeKuttaOrder>::advect(StaggeredGridBasedVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	StaggeredGridBasedVectorField<Dim> forwardField(field.staggeredGrid());
	StaggeredGridBasedVectorField<Dim> backwardField(field.staggeredGrid());
	SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(field, forwardField, flow, dt);
	SemiLagrangianAdvector<Dim, RungeKuttaOrder>::advect(forwardField, backwardField, flow, -dt);
	field.parallelForEach([&](const int axis, const VectorDi &coord) {
		field[axis][coord] = forwardField[axis][coord] + (field[axis][coord] - backwardField[axis][coord]) * real(0.5);
	});
}

template <int Dim, int RungeKuttaOrder>
void MacCormackAdvector<Dim, RungeKuttaOrder>::advect(Particles<Dim> &particles, const VectorField<Dim> &flow, const real dt)
{
	particles.parallelForEach([&](const int i) {
		VectorDr &pos = particles.positions[i];
		const VectorDr forwardPos = trace(pos, flow, dt);
		const VectorDr backwardPos = trace(forwardPos, flow, -dt);
		pos = forwardPos + (pos - backwardPos) * real(0.5);
	});
}

template class EulerianAdvector<2>;
template class EulerianAdvector<3>;

template class SemiLagrangianAdvector<2, 1>;
template class SemiLagrangianAdvector<2, 2>;
template class SemiLagrangianAdvector<2, 3>;
template class SemiLagrangianAdvector<2, 4>;
template class SemiLagrangianAdvector<3, 1>;
template class SemiLagrangianAdvector<3, 2>;
template class SemiLagrangianAdvector<3, 3>;
template class SemiLagrangianAdvector<3, 4>;

template class MacCormackAdvector<2, 1>;
template class MacCormackAdvector<2, 2>;
template class MacCormackAdvector<2, 3>;
template class MacCormackAdvector<2, 4>;
template class MacCormackAdvector<3, 1>;
template class MacCormackAdvector<3, 2>;
template class MacCormackAdvector<3, 3>;
template class MacCormackAdvector<3, 4>;

}
