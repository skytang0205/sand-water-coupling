#include "LevelSetReinitializer.h"

#include <algorithm>
#include <limits>

#include <cmath>

namespace PhysX {

template <int Dim>
FastMarchingReinitializer<Dim>::FastMarchingReinitializer(const Grid<Dim> *const grid, const int reinitMaxIters) :
	LevelSetReinitializer<Dim>(reinitMaxIters),
	_tent(grid),
	_valid(grid)
{ }

template <int Dim>
void FastMarchingReinitializer<Dim>::reinitialize(GridBasedImplicitSurface<Dim>&levelSet)
{
	auto &phi = levelSet.signedDistanceField();
	const real bandWidth = _reinitMaxIters * phi.spacing();
	_tent.setConstant(bandWidth > 0 ? bandWidth : std::numeric_limits<real>::infinity());
	_valid.setZero();
	initInterface(phi);
	fastMarching();
	phi.parallelForEach([&](const VectorDi &cell) {
		phi[cell] = (phi[cell] > 0 ? 1 : -1) * _tent[cell];
	});
}

template <int Dim>
void FastMarchingReinitializer<Dim>::initInterface(const GridBasedScalarField<Dim> &phi)
{
	phi.forEach([&](const VectorDi &cell) {
		for (int axis = 0; axis < Dim; axis++) {
			const VectorDi &nbCell = cell + VectorDi::Unit(axis);
			if (phi.isValid(nbCell) && phi[cell] != phi[nbCell]) {
				if (const auto idx = phi.index(cell); !_valid[cell]) {
					_valid[cell] = true;
				}
				if (const auto idx = phi.index(nbCell); !_valid[nbCell]) {
					_valid[nbCell] = true;
				}
			}
		}
	});
}

template class LevelSetReinitializer<2>;
template class LevelSetReinitializer<3>;

template class FastMarchingReinitializer<2>;
template class FastMarchingReinitializer<3>;

}
