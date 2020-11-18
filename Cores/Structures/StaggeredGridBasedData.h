#pragma once

#include "Structures/GridBasedData.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/IO.h"

namespace PhysX {

template <int Dim, typename Type = real>
class StaggeredGridBasedData
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> *_grid = nullptr;
	std::array<GridBasedData<Dim>, Dim> _components;

public:

	StaggeredGridBasedData(const StaggeredGrid<Dim> *const grid, const Type &value = Type()) { resize(grid, value); }

	StaggeredGridBasedData() = default;
	virtual ~StaggeredGridBasedData() = default;

	void resize(const StaggeredGrid<Dim> *const grid, const Type &value = Type())
	{
		_grid = grid;
		for (int axis = 0; axis < Dim; axis++)
			_components[axis].resize(_grid->faceGrid(axis), value);
	}

	bool isBoundary(const int axis, const VectorDi &face) { return _grid->isBoundaryFace(axis, face); }
	real spacing() const { return _grid->spacing(); }

	GridBasedData<Dim> &operator[](const int axis) { return _components[axis]; }
	const GridBasedData<Dim> &operator[](const int axis) const { return _components[axis]; }

	Type min() const
	{
		if constexpr (Dim == 2) return std::min(_components[0].min(), _components[1].min());
		else return std::min({ _components[0].min(), _components[1].min(), _components[2].min() });
	}

	Type max() const
	{
		if constexpr (Dim == 2) return std::max(_components[0].max(), _components[1].max());
		else return std::max({ _components[0].max(), _components[1].max(), _components[2].max() });
	}

	Type absoluteMax() const
	{
		if constexpr (Dim == 2) return std::max(_components[0].absoluteMax(), _components[1].absoluteMax());
		else return std::max({ _components[0].absoluteMax(), _components[1].absoluteMax(), _components[2].absoluteMax() });
	}

	void forEach(const std::function<void(const int, const VectorDi &)> &func) const { _grid->forEachFace(func); }
	void parallelForEach(const std::function<void(const int, const VectorDi &)> &func) const { _grid->parallelForEachFace(func); }

	void load(std::istream &in) { for (int axis = 0; axis < Dim; axis++) _components[axis].load(in); }
	void save(std::ostream &out) const { for (int axis = 0; axis < Dim; axis++) _components[axis].save(out); }
};

}
