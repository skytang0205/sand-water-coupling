#include "Grid.h"

namespace PhysX {

template <int Dim>
auto Grid<Dim>::linearNearbyDataPoints(const VectorDr &pos) const->std::array<VectorDi, _kCntNb1>
{
	const VectorDi lower = getLinearLower(pos);

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
auto Grid<Dim>::cubicNearbyDataPoints(const VectorDr &pos) const->std::array<VectorDi, _kCntNb3>
{
	const VectorDi lower = getCubicLower(pos);

	std::array<VectorDi, _kCntNb3> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				dataPoints[j << 2 | i] = lower + Vector2i(i, j);
	}
	else {
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					dataPoints[k << 4 | j << 2 | i] = lower + Vector3i(i, j, k);
	}
	return dataPoints;
}

template <int Dim>
auto Grid<Dim>::linearIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, _kCntNb1>
{
	const VectorDi lower = getLinearLower(pos);
	const VectorDr frac = getLowerFrac(pos, lower);
	const std::array<VectorDr, 2> w = {
		VectorDr::Ones() - frac,
		frac
	};

	if constexpr (Dim == 2) {
		return {
			IntrplDataPoint(lower + Vector2i(0, 0), w[0][0] * w[0][1]),
			IntrplDataPoint(lower + Vector2i(1, 0), w[1][0] * w[0][1]),
			IntrplDataPoint(lower + Vector2i(0, 1), w[0][0] * w[1][1]),
			IntrplDataPoint(lower + Vector2i(1, 1), w[1][0] * w[1][1])
		};
	}
	else {
		return {
			IntrplDataPoint(lower + Vector3i(0, 0, 0), w[0][0] * w[0][1] * w[0][2]),
			IntrplDataPoint(lower + Vector3i(1, 0, 0), w[1][0] * w[0][1] * w[0][2]),
			IntrplDataPoint(lower + Vector3i(0, 1, 0), w[0][0] * w[1][1] * w[0][2]),
			IntrplDataPoint(lower + Vector3i(1, 1, 0), w[1][0] * w[1][1] * w[0][2]),
			IntrplDataPoint(lower + Vector3i(0, 0, 1), w[0][0] * w[0][1] * w[1][2]),
			IntrplDataPoint(lower + Vector3i(1, 0, 1), w[1][0] * w[0][1] * w[1][2]),
			IntrplDataPoint(lower + Vector3i(0, 1, 1), w[0][0] * w[1][1] * w[1][2]),
			IntrplDataPoint(lower + Vector3i(1, 1, 1), w[1][0] * w[1][1] * w[1][2])
		};
	}
}

template <int Dim>
auto Grid<Dim>::gradientLinearIntrplDataPoints(const VectorDr &pos) const->std::array<GradientIntrplDataPoint, _kCntNb1>
{
	const VectorDi lower = getLinearLower(pos);
	const VectorDr frac = getLowerFrac(pos, lower);
	const real invDx = 1 / _spacing;
	const std::array<VectorDr, 2> w = {
		VectorDr::Ones() - frac,
		frac
	};

	if constexpr (Dim == 2) {
		return {
			GradientIntrplDataPoint(lower + Vector2i(0, 0), Vector2r(-w[0][1], -w[0][0]) * invDx),
			GradientIntrplDataPoint(lower + Vector2i(1, 0), Vector2r( w[0][1], -w[1][0]) * invDx),
			GradientIntrplDataPoint(lower + Vector2i(0, 1), Vector2r(-w[1][1],  w[0][0]) * invDx),
			GradientIntrplDataPoint(lower + Vector2i(1, 1), Vector2r( w[1][1],  w[1][0]) * invDx)
		};
	}
	else {
		return {
			GradientIntrplDataPoint(lower + Vector3i(0, 0, 0), Vector3r(-w[0][1] * w[0][2], -w[0][0] * w[0][2], -w[0][0] * w[0][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(1, 0, 0), Vector3r( w[0][1] * w[0][2], -w[1][0] * w[0][2], -w[1][0] * w[0][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(0, 1, 0), Vector3r(-w[1][1] * w[0][2],  w[0][0] * w[0][2], -w[0][0] * w[1][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(1, 1, 0), Vector3r( w[1][1] * w[0][2],  w[1][0] * w[0][2], -w[1][0] * w[1][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(0, 0, 1), Vector3r(-w[0][1] * w[1][2], -w[0][0] * w[1][2],  w[0][0] * w[0][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(1, 0, 1), Vector3r( w[0][1] * w[1][2], -w[1][0] * w[1][2],  w[1][0] * w[0][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(0, 1, 1), Vector3r(-w[1][1] * w[1][2],  w[0][0] * w[1][2],  w[0][0] * w[1][1]) * invDx),
			GradientIntrplDataPoint(lower + Vector3i(1, 1, 1), Vector3r( w[1][1] * w[1][2],  w[1][0] * w[1][2],  w[1][0] * w[1][1]) * invDx)
		};
	}
}

template <int Dim>
auto Grid<Dim>::quadraticBasisSplineIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, _kCntNb2>
{
	const VectorDi lower = getQuadraticLower(pos);
	const VectorDr frac = getLowerFrac(pos, lower);
	const std::array<VectorDr, 3> w = {
		(frac - VectorDr::Ones() * 0).unaryExpr(MathFunc::cubicBasisSplineCoefficient),
		(frac - VectorDr::Ones() * 1).unaryExpr(MathFunc::cubicBasisSplineCoefficient),
		(frac - VectorDr::Ones() * 2).unaryExpr(MathFunc::cubicBasisSplineCoefficient)
	};

	std::array<IntrplDataPoint, _kCntNb2> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 3; j++)
			for (int i = 0; i < 3; i++)
				dataPoints[(j << 1) + j + i] = IntrplDataPoint(
					lower + Vector2i(i, j),
					w[i][0] * w[j][1]
				);
	}
	else {
		for (int k = 0; k < 3; k++)
			for (int j = 0; j < 3; j++)
				for (int i = 0; i < 3; i++)
					dataPoints[(k << 3) + k + (j << 1) + j + i] = IntrplDataPoint(
						lower + Vector3i(i, j, k),
						w[i][0] * w[j][1] * w[k][2]
					);
	}
	return dataPoints;
}

template <int Dim>
auto Grid<Dim>::cubicBasisSplineIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, _kCntNb3>
{
	const VectorDi lower = getCubicLower(pos);
	const VectorDr frac = getLowerFrac(pos, lower);
	const std::array<VectorDr, 4> w = {
		(frac - VectorDr::Ones() * 0).unaryExpr(MathFunc::cubicBasisSplineCoefficient),
		(frac - VectorDr::Ones() * 1).unaryExpr(MathFunc::cubicBasisSplineCoefficient),
		(frac - VectorDr::Ones() * 2).unaryExpr(MathFunc::cubicBasisSplineCoefficient),
		(frac - VectorDr::Ones() * 3).unaryExpr(MathFunc::cubicBasisSplineCoefficient)
	};

	std::array<IntrplDataPoint, _kCntNb3> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				dataPoints[j << 2 | i] = IntrplDataPoint(
					lower + Vector2i(i, j),
					w[i][0] * w[j][1]
				);
	}
	else {
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					dataPoints[k << 4 | j << 2 | i] = IntrplDataPoint(
						lower + Vector3i(i - 1, j - 1, k - 1),
						w[i][0] * w[j][1] * w[k][2]
					);
	}
	return dataPoints;
}

template <int Dim>
auto Grid<Dim>::cubicCatmullRomIntrplDataPoints(const VectorDr &pos) const->std::array<IntrplDataPoint, _kCntNb3>
{
	const VectorDi lower = getCubicLower(pos);
	const VectorDr frac = getLowerFrac(pos, lower);
	const std::array<VectorDr, 4> w = {
		frac.unaryExpr(std::bind(MathFunc::cubicCatmullRomCoefficient, 0, std::placeholders::_1)),
		frac.unaryExpr(std::bind(MathFunc::cubicCatmullRomCoefficient, 1, std::placeholders::_1)),
		frac.unaryExpr(std::bind(MathFunc::cubicCatmullRomCoefficient, 2, std::placeholders::_1)),
		frac.unaryExpr(std::bind(MathFunc::cubicCatmullRomCoefficient, 3, std::placeholders::_1))
	};

	std::array<IntrplDataPoint, _kCntNb3> dataPoints;
	if constexpr (Dim == 2) {
		for (int j = 0; j < 4; j++)
			for (int i = 0; i < 4; i++)
				dataPoints[j << 2 | i] = IntrplDataPoint(
					lower + Vector2i(i, j),
					w[i][0] * w[j][1]
				);
	}
	else {
		for (int k = 0; k < 4; k++)
			for (int j = 0; j < 4; j++)
				for (int i = 0; i < 4; i++)
					dataPoints[k << 4 | j << 2 | i] = IntrplDataPoint(
						lower + Vector3i(i, j, k),
						w[i][0] * w[j][1] * w[k][2]
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
