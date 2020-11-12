#include "Grid.h"

namespace PhysX {

template <int Dim>
void Grid<Dim>::getLerp(const VectorDr &dataOrigin, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const
{
	VectorDi lower = ((pos - dataOrigin) / _spacing).cast<int>();
	VectorDr frac = (pos - dataOrigin - lower.cast<real>() * _spacing) / _spacing;
	if constexpr (Dim == 2) {
		coords = { lower, lower + Vector2i(0, 1), lower + Vector2i(1, 0), lower + Vector2i(1, 1) };
		weights = { (1 - frac.x()) * (1 - frac.y()), (1 - frac.x()) * frac.y(), frac.x() * (1 - frac.y()), frac.x() * frac.y() };
	}
	else {
		coords = {
			lower, lower + Vector3i(0, 0, 1), lower + Vector3i(0, 1, 0), lower + Vector3i(0, 1, 1),
			lower + Vector3i(1, 0, 0), lower + Vector3i(1, 0, 1), lower + Vector3i(1, 1, 0), lower + Vector3i(1, 1, 1)
		};
		weights = {
			(1 - frac.x()) * (1 - frac.y()) * (1 - frac.z()), (1 - frac.x()) * (1 - frac.y()) * frac.z(),
			(1 - frac.x()) * frac.y() * (1 - frac.z()), (1 - frac.x()) * frac.y() * frac.z(),
			frac.x() * (1 - frac.y()) * (1 - frac.z()), frac.x() * (1 - frac.y()) * frac.z(),
			frac.x() * frac.y() * (1 - frac.z()), frac.x() * frac.y() * frac.z()
		};
	}
}

template class Grid<2>;
template class Grid<3>;

}
