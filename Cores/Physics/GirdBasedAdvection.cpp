#include "GirdBasedAdvection.h"

namespace PhysX {

template<int Dim>
void GridBasedAdvection<Dim>::advect(VectorGridField<Dim, FaceCentered> &field, const VectorField<Dim> &flow, const real dt)
{
	auto newField = field;
}

template class GridBasedAdvection<2>;
template class GridBasedAdvection<3>;

}
