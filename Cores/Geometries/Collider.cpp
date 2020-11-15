#include "Collider.h"

namespace PhysX {

template class Collider<2>;
template class Collider<3>;

template class StaticCollider<2>;
template class StaticCollider<3>;

template class DynamicCollider<2>;
template class DynamicCollider<3>;

}
