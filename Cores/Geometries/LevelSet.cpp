#include "LevelSet.h"

namespace PhysX {

template <int Dim>
void LevelSet<Dim>::unionSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::min(_signedDistanceField[coord], surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
void LevelSet<Dim>::intersectSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::max(_signedDistanceField[coord], surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
void LevelSet<Dim>::exceptSurface(const Surface<Dim> &surface)
{
	_signedDistanceField.parallelForEach([&](const VectorDi &coord) {
		_signedDistanceField[coord] = std::max(_signedDistanceField[coord], -surface.signedDistance(_signedDistanceField.position(coord)));
	});
}

template <int Dim>
real LevelSet<Dim>::curvature(const VectorDr &pos) const
{
	const real dx = _signedDistanceField.spacing();
	const real invDx = _signedDistanceField.invSpacing();
	real acc = 0;
	for (int i = 0; i < Dim; i++) {
		acc += closestNormal(pos + VectorDr::Unit(i) * dx / 2)[i] - closestNormal(pos - VectorDr::Unit(i) * dx / 2)[i];
	}
	acc *= invDx;
	return std::abs(acc) < invDx ? acc : (acc < 0 ? -1 : 1) * invDx;
}

template class LevelSet<2>;
template class LevelSet<3>;

}
