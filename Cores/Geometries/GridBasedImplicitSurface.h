#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Structures/GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
class GridBasedImplicitSurface final : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	GridBasedScalarField<Dim> _signedDistanceField;

public:

	GridBasedImplicitSurface(const Grid<Dim> *const grid) : _signedDistanceField(grid) { }
	GridBasedImplicitSurface(const Grid<Dim> *const grid, const Surface<Dim> &surface) : _signedDistanceField(grid) { setFromSurface(surface); }

	virtual ~GridBasedImplicitSurface() = default;

	void setFromSurface(const Surface<Dim> &surface);

	GridBasedScalarField<Dim> &signedDistanceField() { return _signedDistanceField; }
	const GridBasedScalarField<Dim> &signedDistanceField() const { return _signedDistanceField; }

	virtual VectorDr closestNormal(const VectorDr &pos) const
	{
		VectorDr n = _signedDistanceField.gradient(pos);
		if (n.any()) n.normalize();
		return n;
	}

	virtual real signedDistance(const VectorDr &pos) const { return _signedDistanceField(pos); }
};

}
