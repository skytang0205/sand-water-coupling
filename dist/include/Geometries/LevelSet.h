#pragma once

#include "Geometries/ImplicitSurface.h"
#include "Structures/GridBasedScalarField.h"

namespace PhysX {

template <int Dim>
class LevelSet : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	GridBasedScalarField<Dim> _signedDistanceField;

public:

	LevelSet(const Grid<Dim> *const grid) : _signedDistanceField(grid, std::numeric_limits<real>::infinity()) { }
	LevelSet(const Grid<Dim> *const grid, const Surface<Dim> &surface) : _signedDistanceField(grid, std::numeric_limits<real>::infinity()) { unionSurface(surface); }

	virtual ~LevelSet() = default;

	void clear() { _signedDistanceField.setConstant(std::numeric_limits<real>::infinity()); }
	void unionSurface(const Surface<Dim> &surface);
	void intersectSurface(const Surface<Dim> &surface);
	void exceptSurface(const Surface<Dim> &surface);

	GridBasedScalarField<Dim> &signedDistanceField() { return _signedDistanceField; }
	const GridBasedScalarField<Dim> &signedDistanceField() const { return _signedDistanceField; }

	virtual VectorDr closestNormal(const VectorDr &pos) const { return (_signedDistanceField.gradient(pos)).normalized(); }
	virtual real signedDistance(const VectorDr &pos) const { return _signedDistanceField(pos); }

	real curvature(const VectorDr &pos) const;
};

}
