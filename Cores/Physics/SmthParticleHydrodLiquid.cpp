#include "SmthParticleHydrodLiquid.h"

namespace PhysX {

template <int Dim>
void SmthParticleHydrodLiquid<Dim>::writeDescription(YAML::Node &root) const
{
}

template class SmthParticleHydrodLiquid<2>;
template class SmthParticleHydrodLiquid<3>;

}
