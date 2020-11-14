#pragma once

#include "Structures/Grid.h"
#include "Utilities/IO.h"

#include <vector>

namespace PhysX {

template <int Dim, typename Type = real>
class GridBasedData
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const Grid<Dim> *_grid = nullptr;
	std::vector<Type> _data;

public:

	GridBasedData(const Grid<Dim> *const grid, const Type &value = Type()) { resize(grid, value); }

	GridBasedData() = default;
	virtual ~GridBasedData() = default;

	void resize(const Grid<Dim> *const grid, const Type &value = Type())
	{
		_grid = grid;
		_data.resize(_grid->dataCount(), value);
	}

	VectorDi size() const { return _grid->dataSize(); }
	size_t count() const { return _data.size(); }
	VectorDr position(const VectorDi &coord) const { return _grid->dataPosition(coord); }

	const Type &operator[](const VectorDi &coord) const { return _data[_grid->index(coord)]; }
	Type &operator[](const VectorDi &coord) { return _data[_grid->index(coord)]; }

	void forEach(const std::function<void(const VectorDi &)> &func) const { _grid->forEach(func); }
	void parallelForEach(const std::function<void(const VectorDi &)> &func) const { _grid->parallelForEach(func); }

	void read(std::istream &in) { IO::readArray(in, _data.data(), _data.size()); }
	void write(std::ostream &out) const { IO::writeArray(out, _data.data(), _data.size()); }
};

}
