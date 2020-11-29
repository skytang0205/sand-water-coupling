#include "Grid.h"

#include "Utilities/MathFunc.h"

namespace PhysX {

template <int Dim>
auto Grid<Dim>::oneLayerNearbyDataPoints(const VectorDr &pos) const->std::array<VectorDi, 1 << Dim>
{
	const VectorDi lower = getLowerCoord(pos);
	if constexpr (Dim == 2) {
		return {
			lower + Vector2i(0, 0), lower + Vector2i(1, 0),
			lower + Vector2i(0, 1), lower + Vector2i(1, 1)
		};
	}
	else {
		return {
			lower + Vector3i(0, 0, 0), lower + Vector3i(1, 0, 0),
			lower + Vector3i(0, 1, 0), lower + Vector3i(1, 1, 0),
			lower + Vector3i(0, 0, 1), lower + Vector3i(1, 0, 1),
			lower + Vector3i(0, 1, 1), lower + Vector3i(1, 1, 1)
		};
	}
}

template <int Dim>
auto Grid<Dim>::twoLayersNearbyDataPoints(const VectorDr &pos) const->std::array<VectorDi, 1 << (Dim << 1)>
{
	const VectorDi lower = getLowerCoord(pos);
	std::array<VectorDi, 1 << (Dim << 1)> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				dataPoints[j << 2 | i] = lower + Vector2i(i - 1, j - 1);
	}
	else {
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					dataPoints[k << 4 | j << 2 | i] = lower + Vector3i(i - 1, j - 1, k - 1);
	}
	return dataPoints;
}

template <int Dim>
auto Grid<Dim>::linearIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, 1 << Dim>
{
	const auto [lower, frac] = getLowerCoordAndFrac(pos);
	if constexpr (Dim == 2) {
		return {
			IntrplDataPoint(lower + Vector2i(0, 0), (1 - frac.x()) * (1 - frac.y())),
			IntrplDataPoint(lower + Vector2i(1, 0), frac.x() * (1 - frac.y())),
			IntrplDataPoint(lower + Vector2i(0, 1), (1 - frac.x()) * frac.y()),
			IntrplDataPoint(lower + Vector2i(1, 1), frac.x() * frac.y())
		};
	}
	else {
		return {
			IntrplDataPoint(lower + Vector3i(0, 0, 0), (1 - frac.x()) * (1 - frac.y()) * (1 - frac.z())),
			IntrplDataPoint(lower + Vector3i(1, 0, 0), frac.x() * (1 - frac.y()) * (1 - frac.z())),
			IntrplDataPoint(lower + Vector3i(0, 1, 0), (1 - frac.x()) * frac.y() * (1 - frac.z())),
			IntrplDataPoint(lower + Vector3i(1, 1, 0), frac.x() * frac.y() * (1 - frac.z())),
			IntrplDataPoint(lower + Vector3i(0, 0, 1), (1 - frac.x()) * (1 - frac.y()) * frac.z()),
			IntrplDataPoint(lower + Vector3i(1, 0, 1), frac.x() * (1 - frac.y()) * frac.z()),
			IntrplDataPoint(lower + Vector3i(0, 1, 1), (1 - frac.x()) * frac.y() * frac.z()),
			IntrplDataPoint(lower + Vector3i(1, 1, 1), frac.x() * frac.y() * frac.z())
		};
	}
}

template <int Dim>
auto Grid<Dim>::cubicCatmullRomIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, 1 << (Dim << 1)>
{
	const auto [lower, frac] = getLowerCoordAndFrac(pos);
	std::array<IntrplDataPoint, 1 << (Dim << 1)> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				dataPoints[j << 2 | i] = IntrplDataPoint(
					lower + Vector2i(i - 1, j - 1),
					MathFunc::cubicCatmullRomSplineCoefficient(i, frac[0]) * MathFunc::cubicCatmullRomSplineCoefficient(j, frac[1])
				);
	}
	else {
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					dataPoints[k << 4 | j << 2 | i] = IntrplDataPoint(
						lower + Vector3i(i - 1, j - 1, k - 1),
						MathFunc::cubicCatmullRomSplineCoefficient(i, frac[0]) * MathFunc::cubicCatmullRomSplineCoefficient(j, frac[1]) * MathFunc::cubicCatmullRomSplineCoefficient(k, frac[2])
					);
	}
	return dataPoints;
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
