#pragma once

#include "Surface.h"

template <int Dim>
class SurfaceMesh : public Surface<Dim>
{
	static_assert(Dim >= 2 && Dim <= 3, "Dimension can only be 2 or 3.");

	DECLARE_DIM_TYPES(Dim)

protected:


public:

};
