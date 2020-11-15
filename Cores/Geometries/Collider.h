#pragma once

#include "Geometries/Surface.h"

#include <memory>

namespace PhysX {

template <int Dim>
class Collider
{
	DECLARE_DIM_TYPES(Dim)

protected:

	real _restitutionCoefficient;
	real _frictionCoefficient;

public:

	Collider(const real restitutionCoefficient = 1, const real frictionCoefficient = 0) :
		_restitutionCoefficient(restitutionCoefficient),
		_frictionCoefficient(frictionCoefficient)
	{ }

	Collider(const Collider &rhs) = delete;
	Collider &operator=(const Collider &rhs) = delete;
	virtual ~Collider() = default;

	virtual Surface<Dim> *surface() const = 0;

	virtual VectorDr velocityAt(const VectorDr &pos) const = 0;
	real restitutionCoefficient() const { return _restitutionCoefficient; }
	real frictionCoefficient() const { return _frictionCoefficient; }
};

template <int Dim>
class StaticCollider : public Collider<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const std::unique_ptr<Surface<Dim>> _surface;

public:

	StaticCollider(std::unique_ptr<Surface<Dim>> surface, const real restitutionCoefficient = 1, const real frictionCoefficient = 0) :
		Collider<Dim>(restitutionCoefficient, frictionCoefficient),
		_surface(std::move(surface))
	{ }

	StaticCollider(const StaticCollider &rhs) = delete;
	StaticCollider &operator=(const StaticCollider &rhs) = delete;
	virtual ~StaticCollider() = default;

	virtual VectorDr velocityAt(const VectorDr &pos) const override { return VectorDr::Zero(); }
	virtual Surface<Dim> *surface() const override final { return _surface.get(); }
};

template <int Dim>
class DynamicCollider : public Collider<Dim>
{
};

}
