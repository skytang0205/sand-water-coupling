#include "Grid.h"

#include "Utilities/MathFunc.h"

namespace PhysX {

template <int Dim>
void Grid<Dim>::getLerpCoordsAndWeights(const VectorDr &pos, std::array<VectorDi, 1 << Dim> &coords, std::array<real, 1 << Dim> &weights) const
{
	VectorDi lower;
	VectorDr frac;
	getLowerCoordAndFrac(pos, lower, frac);
	if constexpr (Dim == 2) {
		coords = { lower, lower + Vector2i(1, 0), lower + Vector2i(0, 1), lower + Vector2i(1, 1) };
		weights = { (1 - frac.x()) * (1 - frac.y()), frac.x() * (1 - frac.y()), (1 - frac.x()) * frac.y(), frac.x() * frac.y() };
	}
	else {
		coords = {
			lower, lower + Vector3i(1, 0, 0), lower + Vector3i(0, 1, 0), lower + Vector3i(1, 1, 0),
			lower + Vector3i(0, 0, 1), lower + Vector3i(1, 0, 1), lower + Vector3i(0, 1, 1), lower + Vector3i(1, 1, 1)
		};
		weights = {
			(1 - frac.x()) * (1 - frac.y()) * (1 - frac.z()), frac.x() * (1 - frac.y()) * (1 - frac.z()),
			(1 - frac.x()) * frac.y() * (1 - frac.z()), frac.x() * frac.y() * (1 - frac.z()),
			(1 - frac.x()) * (1 - frac.y()) * frac.z(), frac.x() * (1 - frac.y()) * frac.z(),
			(1 - frac.x()) * frac.y() * frac.z(), frac.x() * frac.y() * frac.z()
		};
	}
}

template <int Dim>
void Grid<Dim>::getSplineCoordsAndWeights(const VectorDr &pos, std::array<VectorDi, 1 << (Dim << 1)> &coords, std::array<real, 1 << (Dim << 1)> &weights) const
{
	VectorDi lower;
	VectorDr frac;
	getLowerCoordAndFrac(pos, lower, frac);
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++) {
			for (int i = 0; i < 4; i++) {
				coords[j << 2 | i] = (lower + Vector2i(i - 1, j - 1)).cwiseMax(Vector2i::Zero()).cwiseMin(_dataSize - Vector2i::Ones());
				weights[j << 2 | i] = MathFunc::cubicSplineCoefficient(i, frac[0]) * MathFunc::cubicSplineCoefficient(j, frac[1]);
			}
		}
	}
	else {
		for (int k = 0; k < 4; k++) {
			for (int j = 0; j < 4; j++) {
				for (int i = 0; i < 4; i++) {
					coords[k << 4 | j << 2 | i] = (lower + Vector3i(i - 1, j - 1, k - 1)).cwiseMax(Vector3i::Zero()).cwiseMin(_dataSize - Vector3i::Ones());
					weights[k << 4 | j << 2 | i] = MathFunc::cubicSplineCoefficient(i, frac[0]) * MathFunc::cubicSplineCoefficient(j, frac[1]) * MathFunc::cubicSplineCoefficient(k, frac[2]);
				}
			}
		}
	}
}

template <int Dim>
void Grid<Dim>::forEach(const std::function<void(const VectorDi &)> &func) const
{
	if constexpr (Dim == 2) {
		for (int j = 0; j < _dataSize.y(); j++)
			for (int i = 0; i < _dataSize.x(); i++)
				func(VectorDi(i, j));
	}
	else {
		for (int k = 0; k < _dataSize.z(); k++)
			for (int j = 0; j < _dataSize.y(); j++)
				for (int i = 0; i < _dataSize.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template <int Dim>
void Grid<Dim>::parallelForEach(const std::function<void(const VectorDi &)> &func) const
{
	if constexpr (Dim == 2) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int j = 0; j < _dataSize.y(); j++)
			for (int i = 0; i < _dataSize.x(); i++)
				func(VectorDi(i, j));
	}
	else {
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (int k = 0; k < _dataSize.z(); k++)
			for (int j = 0; j < _dataSize.y(); j++)
				for (int i = 0; i < _dataSize.x(); i++)
					func(VectorDi(i, j, k));
	}
}

template class Grid<2>;
template class Grid<3>;

}
