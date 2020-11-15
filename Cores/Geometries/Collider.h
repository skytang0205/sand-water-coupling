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
class StaticCollider final : public Collider<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::unique_ptr<Surface<Dim>> _surface;

public:

	StaticCollider(const real restitutionCoefficient = 1, const real frictionCoefficient = 0, std::unique_ptr<Surface<Dim>> surface) :
		Collider<Dim>(restitutionCoefficient, frictionCoefficient),
		_surface(std::move(surface))
	{ }

	StaticCollider(const StaticCollider &rhs) = delete;
	StaticCollider &operator=(const StaticCollider &rhs) = delete;
	virtual ~StaticCollider() = default;

	virtual Surface<Dim> *surface() const override { return _surface.get(); }
};

}
