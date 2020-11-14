#pragma once

#include "Geometries/Surface.h"

#include <memory>

namespace PhysX {

template <int Dim>
class Collider
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::shared_ptr<Surface<Dim>> _surface;

public:

	Surface<Dim> *surface() const { return _surface.get(); }
};

}
