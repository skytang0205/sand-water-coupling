#include "GirdBasedAdvection.h"

namespace PhysX {

template<int Dim>
void SemiLagrangianAdvection<Dim>::advect(GridBasedStaggeredVectorField<Dim> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
	newField.parallelForEach([&](const int axis, const VectorDi &face) {
			const VectorDr startPos = newField[axis].position(face);
			const VectorDr midPos = startPos - flow(startPos) * dt * real(0.5);
			newField[axis][face] = field[axis](startPos - flow(midPos) * dt);
		});
	field = newField;
}

template class GridBasedAdvection<2>;
template class GridBasedAdvection<3>;

template class SemiLagrangianAdvection<2>;
template class SemiLagrangianAdvection<3>;

}
