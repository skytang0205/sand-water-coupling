#pragma once

#include "Geometries/Collider.h"
#include "Structures/StaggeredGridBasedVectorField.h"

namespace PhysX {

template <int Dim>
class EulerianProjector
{
public:

	EulerianProjector() = default;
	EulerianProjector(const EulerianProjector &rhs) = delete;
	EulerianProjector &operator=(const EulerianProjector &rhs) = delete;
	virtual ~EulerianProjector() = default;

	void project(StaggeredGridBasedVectorField<Dim> &velocity, const Surface<Dim> *const fluidSurface, const Collider<Dim> *const collider);

protected:

};

}
