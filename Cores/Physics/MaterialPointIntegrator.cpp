#include "MaterialPointIntegrator.h"

namespace PhysX {

template <int Dim>
void MpSemiImplicitIntegrator<Dim>::integrate(
	GridBasedVectorData<Dim> &velocity,
	const GridBasedScalarData<Dim> &mass,
	const std::vector<std::unique_ptr<MaterialPointSubstance<Dim>>> &substances,
	const real dt,
	const GridBasedData<Dim, uchar> &collided)
{
	// TODO: Impl
}

template class MaterialPointIntegrator<2>;
template class MaterialPointIntegrator<3>;
template class MpSymplecticEulerIntegrator<2>;
template class MpSymplecticEulerIntegrator<3>;
template class MpSemiImplicitIntegrator<2>;
template class MpSemiImplicitIntegrator<3>;

}
