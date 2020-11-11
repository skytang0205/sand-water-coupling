#pragma once

#include "Structures/Field.h"

namespace PhysX {

template <int Dim>
class Grid
{
	static_assert(2 <= Dim && Dim <= 3, "Dimension must be 2 or 3.");
	DECLARE_DIM_TYPES(Dim)

protected:

	const real _spacing;
	const VectorDi _resolution;
	const VectorDr _origin;

public:

	Grid(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero()) :
		_spacing(spacing),
		_resolution(resolution),
		_origin(center - resolution.cast<real>() * spacing * real(0.5))
	{ }

	Grid(const Grid &rhs) = default;
	Grid &operator=(const Grid &rhs) = default;
	virtual ~Grid() = default;

	VectorDr cellCenter(const VectorDi &coord) { return _origin + (coord.cast<real>() + VectorDr::Ones() * real(0.5)) * _spacing; }
	VectorDr faceCenter(const int axis, const VectorDi &coord) { return _origin + (coord.cast<real>() + (VectorDr::Ones() - VectorDr::Unit(axis)) * real(0.5)) * _spacing; }
};

template <int Dim>
class ScalarGridField : public Grid<Dim>, public ScalarField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	ScalarGridField(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero()) :
		Grid<Dim>(spacing, resolution, center)
	{ }

	ScalarGridField() = delete;
	ScalarGridField(const ScalarGridField &rhs) = default;
	ScalarGridField &operator=(const ScalarGridField &rhs) = default;
	virtual ~ScalarGridField() = default;

	virtual const real &operator[](const VectorDi &coord) const = 0;
	virtual real &operator[](const VectorDi &coord) = 0;
};

template <int Dim>
class VectorGridField : public Grid<Dim>, public VectorField<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	VectorGridField(const real spacing, const VectorDi &resolution, const VectorDr &center = VectorDr::Zero()) :
		Grid<Dim>(spacing, resolution, center)
	{ }

	VectorGridField() = delete;
	VectorGridField(const VectorGridField &rhs) = default;
	VectorGridField &operator=(const VectorGridField &rhs) = default;
	virtual ~VectorGridField() = default;
};

}
