#pragma once

#include "Utilities/MathFunc.h"
#include "Utilities/Types.h"

#include <array>
#include <functional>

namespace PhysX {

template <int Dim>
class Grid final
{
	DECLARE_DIM_TYPES(Dim)

	using IntrplDataPoint = std::pair<VectorDi, real>;
	using GradientIntrplDataPoint = std::pair<VectorDi, VectorDr>;

protected:

	static constexpr int _kCntNb1 = MathFunc::pow(2, Dim);
	static constexpr int _kCntNb2 = MathFunc::pow(3, Dim);
	static constexpr int _kCntNb3 = MathFunc::pow(4, Dim);

	const real _spacing;
	const VectorDi _dataSize;
	const VectorDr _dataOrigin;

public:

	Grid(const real spacing, const VectorDi &dataSize, const VectorDr &dataOrigin) :
		_spacing(spacing),
		_dataSize(dataSize),
		_dataOrigin(dataOrigin)
	{ }

	Grid &operator=(const Grid &rhs) = delete;
	virtual ~Grid() = default;

	bool isInside(const VectorDi &coord, const int offset) const { return (coord.array() >= offset).all() && (coord.array() < _dataSize.array() - offset).all(); }
	bool isValid(const VectorDi &coord) const { return isInside(coord, 0); }

	real spacing() const { return _spacing; }
	VectorDi dataSize() const { return _dataSize; }
	VectorDr dataOrigin() const { return _dataOrigin; }

	size_t dataCount() const { return _dataSize.template cast<size_t>().prod(); }
	VectorDr dataPosition(const VectorDi &coord) const { return _dataOrigin + coord.template cast<real>() * _spacing; }
	VectorDi clamp(const VectorDi &coord) const { return coord.cwiseMax(0).cwiseMin(_dataSize - VectorDi::Ones()); }

	size_t index(const VectorDi &coord) const
	{
		if constexpr (Dim == 2) return coord.x() + size_t(_dataSize.x()) * coord.y();
		else return coord.x() + _dataSize.x() * (coord.y() + size_t(_dataSize.y()) * coord.z());
	}

	VectorDi coordinate(const size_t index) const
	{
		if constexpr (Dim == 2) return VectorDi(index % _dataSize.x(), index / _dataSize.x());
		else return VectorDi(
			int(index % _dataSize.x()),
			int(index / _dataSize.x() % _dataSize.y()),
			int(index / _dataSize.x() / _dataSize.y()));
	}

	VectorDi getLinearLower(const VectorDr &pos) const { return ((pos - _dataOrigin) / _spacing).array().floor().template cast<int>().matrix(); }
	VectorDi getQuadraticLower(const VectorDr &pos) const { return ((pos - _dataOrigin) / _spacing - VectorDr::Ones() / 2).array().floor().template cast<int>().matrix(); }
	VectorDi getCubicLower(const VectorDr &pos) const { return ((pos - _dataOrigin) / _spacing).array().floor().template cast<int>().matrix() - VectorDi::Ones(); }
	VectorDr getLowerFrac(const VectorDr &pos, const VectorDi &lower) const { return (pos - _dataOrigin - lower.template cast<real>() * _spacing) / _spacing; }

	std::array<VectorDi, _kCntNb1> linearNearbyDataPoints(const VectorDr &pos) const;
	std::array<VectorDi, _kCntNb3> cubicNearbyDataPoints(const VectorDr &pos) const;

	std::array<IntrplDataPoint, _kCntNb1> linearIntrplDataPoints(const VectorDr &pos) const;
	std::array<GradientIntrplDataPoint, _kCntNb1> gradientLinearIntrplDataPoints(const VectorDr &pos) const;
	std::array<IntrplDataPoint, _kCntNb2> quadraticBasisSplineIntrplDataPoints(const VectorDr &pos) const;
	std::array<IntrplDataPoint, _kCntNb3> cubicBasisSplineIntrplDataPoints(const VectorDr &pos) const;
	std::array<IntrplDataPoint, _kCntNb3> cubicCatmullRomIntrplDataPoints(const VectorDr &pos) const;

	void forEach(const std::function<void(const VectorDi &)> &func) const;
	void parallelForEach(const std::function<void(const VectorDi &)> &func) const;

	static constexpr int numberOfNeighbors() { return Dim << 1; }
	static constexpr int neighborAxis(const int ord) { return ord >> 1; }
	static constexpr int neighborSide(const int ord) { return ord & 1 ? 1 : -1; }
	static VectorDi neighbor(const VectorDi &coord, const int ord) { return coord + VectorDi::Unit(neighborAxis(ord)) * neighborSide(ord); }
};

}
