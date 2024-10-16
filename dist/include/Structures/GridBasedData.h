#pragma once

#include "Structures/Grid.h"
#include "Utilities/IO.h"

#include <algorithm>
#include <numeric>
#include <vector>

namespace PhysX {

template <int Dim, typename Type>
class GridBasedData
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Grid<Dim> *_grid = nullptr;
	std::vector<Type> _data;

public:

	GridBasedData(const Grid<Dim> *const grid, const Type &value = Zero<Type>()) { resize(grid, value); }

	GridBasedData() = default;
	virtual ~GridBasedData() = default;

	void resize(const Grid<Dim> *const grid, const Type &value = Zero<Type>())
	{
		_grid = grid;
		_data.resize(_grid->dataCount(), value);
	}

	bool isInside(const VectorDi &coord, const int offset) const { return _grid->isInside(coord, offset); }
	bool isValid(const VectorDi &coord) const { return _grid->isValid(coord); }

	Type *data() { return _data.data(); }
	const Type *data() const { return _data.data(); }

	const Grid<Dim> *grid() const { return _grid; }
	real spacing() const { return _grid->spacing(); }
	real invSpacing() const { return _grid->invSpacing(); }
	VectorDi size() const { return _grid->dataSize(); }
	VectorDr origin() const { return _grid->dataOrigin(); }
	size_t count() const { return _data.size(); }
	VectorDr position(const VectorDi &coord) const { return _grid->dataPosition(coord); }

	size_t index(const VectorDi &coord) const { return _grid->index(coord); }
	VectorDi coordinate(const size_t index) const { return _grid->coordinate(index); }

	Type &operator[](const size_t index) { return _data[index]; }
	const Type &operator[](const size_t index) const { return _data[index]; }
	Type &operator[](const VectorDi &coord) { return _data[_grid->index(coord)]; }
	const Type &operator[](const VectorDi &coord) const { return _data[_grid->index(coord)]; }
	const Type &at(const VectorDi &coord) const { return _data[_grid->index(_grid->clamp(coord))]; }

	void setConstant(const Type &value) { std::fill(_data.begin(), _data.end(), value); }
	void setZero() { setConstant(Zero<Type>()); }

	template <typename AccType = Type>
	AccType sum() const { return std::accumulate(_data.begin(), _data.end(), Zero<AccType>()); }

	Type min() const { return *std::min_element(_data.begin(), _data.end()); }
	Type max() const { return *std::max_element(_data.begin(), _data.end()); }

	Type absoluteMax() const
	{
		auto minmax = std::minmax_element(_data.begin(), _data.end());
		return std::max(std::abs(*minmax.first), std::abs(*minmax.second));
	}

	real normMax() const
	{
		if constexpr (HasSquaredNorm<Type>) {
			real squaredNormMax = 0;
			for (const auto &val : _data) {
				squaredNormMax = std::max(squaredNormMax, val.squaredNorm());
			}
			return std::sqrt(squaredNormMax);
		}
		else return absoluteMax();
	}

	auto asVectorXr() { return Eigen::Map<VectorXr, Eigen::Aligned>(reinterpret_cast<real *>(_data.data()), _data.size() * (sizeof(Type) / sizeof(real))); }
	auto asVectorXr() const { return Eigen::Map<const VectorXr, Eigen::Aligned>(reinterpret_cast<const real *>(_data.data()), _data.size() * (sizeof(Type) / sizeof(real))); }

	void forEach(const std::function<void(const VectorDi &)> &func) const { _grid->forEach(func); }
	void parallelForEach(const std::function<void(const VectorDi &)> &func) const { _grid->parallelForEach(func); }

	void load(std::istream &in) { IO::readArray(in, _data.data(), _data.size()); }
	void save(std::ostream &out) const { IO::writeArray(out, _data.data(), _data.size()); }
};

template <int Dim> using GridBasedScalarData = GridBasedData<Dim, real>;
template <int Dim> using GridBasedVectorData = GridBasedData<Dim, Vector<Dim, real>>;

}
