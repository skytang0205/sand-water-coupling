#include "Grid.h"

namespace PhysX {

template <int Dim>
void Grid<Dim>::getLerp(const VectorDr &dataOrigin, const VectorDi &dataSize, const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const
{
	VectorDi lower = ((pos - dataOrigin) / _spacing).cast<int>().cwiseMax(0).cwiseMin(dataSize - VectorDi::Ones() * 2);
	VectorDr frac = ((pos - dataOrigin - lower.cast<real>() * _spacing) / _spacing).cwiseMax(0).cwiseMin(1);
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

template <int Dim>
void Grid<Dim>::forEach(const VectorDi &dataSize, const std::function<void(const VectorDi &)> &func) const
{
	if constexpr (Dim == 2) {
		for (int j = 0; j < dataSize.y(); j++)
			for (int i = 0; i < dataSize.x(); i++)
				func(VectorDi(i, j));
	}
	else {
		for (int k = 0; k < dataSize.z(); k++)
			for (int j = 0; j < dataSize.y(); j++)
				for (int i = 0; i < dataSize.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template <int Dim>
void Grid<Dim>::parallelForEach(const VectorDi &dataSize, const std::function<void(const VectorDi &)> &func) const
{
	if constexpr (Dim == 2) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int j = 0; j < dataSize.y(); j++)
			for (int i = 0; i < dataSize.x(); i++)
				func(VectorDi(i, j));
	}
	else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int k = 0; k < dataSize.z(); k++)
			for (int j = 0; j < dataSize.y(); j++)
				for (int i = 0; i < dataSize.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template class Grid<2>;
template class Grid<3>;

}
