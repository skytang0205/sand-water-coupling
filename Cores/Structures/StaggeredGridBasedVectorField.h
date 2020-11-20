#pragma once

#include "Structures/GridBasedScalarField.h"
#include "Structures/StaggeredGrid.h"
#include "Utilities/IO.h"

namespace PhysX {

template <int Dim>
class StaggeredGridBasedVectorField final : public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const StaggeredGrid<Dim> *_grid = nullptr;
	std::array<GridBasedScalarField<Dim>, Dim> _components;

public:

	StaggeredGridBasedVectorField(const StaggeredGrid<Dim> *const grid, const VectorDr &value = VectorDr::Zero()) { resize(grid, value); }

	StaggeredGridBasedVectorField() = default;
	virtual ~StaggeredGridBasedVectorField() = default;

	void resize(const StaggeredGrid<Dim> *const grid, const VectorDr &value = VectorDr::Zero())
	{
		_grid = grid;
		for (int axis = 0; axis < Dim; axis++)
			_components[axis].resize(_grid->faceGrid(axis), value[axis]);
	}

	bool isBoundary(const int axis, const VectorDi &face) { return _grid->isBoundaryFace(axis, face); }
	real spacing() const { return _grid->spacing(); }

	GridBasedScalarField<Dim> &operator[](const int axis) { return _components[axis]; }
	const GridBasedScalarField<Dim> &operator[](const int axis) const { return _components[axis]; }

	virtual VectorDr operator()(const VectorDr &pos) const override;

	real divergenceAtCellCenter(const VectorDi &cell) const;
	virtual real divergence(const VectorDr &pos) const override;

	void setConstant(const VectorDr &value) { for (int axis = 0; axis < Dim; axis++) _components[axis].setConstant(value[axis]); }
	void setZero() { setConstant(VectorDr::Zero()); }

	real min() const
	{
		if constexpr (Dim == 2) return std::min(_components[0].min(), _components[1].min());
		else return std::min({ _components[0].min(), _components[1].min(), _components[2].min() });
	}

	real max() const
	{
		if constexpr (Dim == 2) return std::max(_components[0].max(), _components[1].max());
		else return std::max({ _components[0].max(), _components[1].max(), _components[2].max() });
	}

	real absoluteMax() const
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
