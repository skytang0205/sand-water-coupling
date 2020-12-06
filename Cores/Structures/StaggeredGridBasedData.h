#pragma once

#include "Structures/GridBasedData.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/IO.h"

namespace PhysX {

template <int Dim, typename Type>
class StaggeredGridBasedData
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> *_grid = nullptr;
	std::array<GridBasedData<Dim, Type>, Dim> _components;

public:

	StaggeredGridBasedData(const StaggeredGrid<Dim> *const grid, const Type &value = Zero<Type>()) { resize(grid, value); }

	StaggeredGridBasedData() = default;
	virtual ~StaggeredGridBasedData() = default;

	void resize(const StaggeredGrid<Dim> *const grid, const Type &value = Zero<Type>())
	{
		_grid = grid;
		for (int axis = 0; axis < Dim; axis++)
			_components[axis].resize(_grid->faceGrid(axis), value);
	}

	bool isInside(const int axis, const VectorDi &face) const { return _grid->isInsideFace(axis, face); }
	bool isBoundary(const int axis, const VectorDi &face) const { return _grid->isBoundaryFace(axis, face); }
	const StaggeredGrid<Dim> *staggeredGrid() const { return _grid; }
	real spacing() const { return _grid->spacing(); }

	size_t count() const
	{
		if constexpr (Dim == 2) return _components[0].count() + _components[1].count();
		else return _components[0].count() + _components[1].count() + _components[2].count();
	}

	GridBasedData<Dim, Type> &operator[](const int axis) { return _components[axis]; }
	const GridBasedData<Dim, Type> &operator[](const int axis) const { return _components[axis]; }

	void setConstant(const Type &value) { for (int axis = 0; axis < Dim; axis++) _components[axis].setConstant(value); }
	void setZero() { setConstant(Zero<Type>()); }

	template <typename AccType>
	AccType sum() const
	{
		if constexpr (Dim == 2) return _components[0].template sum<AccType>() + _components[1].template sum<AccType>();
		else return _components[0].template sum<AccType>() + _components[1].template sum<AccType>() + _components[2].template sum<AccType>();
	}

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

template <int Dim> using StaggeredGridBasedScalarData = StaggeredGridBasedData<Dim, real>;
template <int Dim> using StaggeredGridBasedVectorData = StaggeredGridBasedData<Dim, Vector<Dim, real>>;

}
