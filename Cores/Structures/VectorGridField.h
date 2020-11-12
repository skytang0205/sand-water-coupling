#pragma once

#include "Structures/Array.h"
#include "Structures/Field.h"
#include "Structures/Grid.h"

namespace PhysX {

template <int Dim, int DataLayout> class VectorGridField;

template <int Dim>
class VectorGridField<Dim, FaceCentered> : public Grid<Dim>, public VectorField<Dim>
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

	using Grid<Dim>::_spacing;
	using Grid<Dim>::faceSize;
	using Grid<Dim>::getCellLerp;
	using Grid<Dim>::getFaceLerp;

protected:

	std::array<Array<Dim>, Dim> _data;

public:

	VectorGridField(
		const real spacing,
		const VectorDi &resolution,
		const VectorDr &center = VectorDr::Zero(),
		const VectorDr &value = VectorDr::Zero())
		:
		Grid<Dim>(spacing, resolution, center)
	{
		for (int axis = 0; axis < Dim; axis++) {
			_data[axis].resize(faceSize(axis), value[axis]);
		}
	}

	VectorGridField(const VectorGridField &rhs) = default;
	VectorGridField &operator=(const VectorGridField & rhs) = default;
	virtual ~VectorGridField() = default;

	const Array<Dim> &operator[](const int axis) const { return _data[axis]; }
	Array<Dim> &operator[](const int axis) { return _data[axis]; }

	virtual VectorDr operator()(const VectorDr &pos) const override final;

	real divergenceAtCellCenter(const VectorDi &cell) const;
	virtual real divergence(const VectorDr &pos) const override final;
};

}
