#pragma once

#include "Geometries/Surface.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class ImplicitSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using Surface<Dim>::isInside;

	ImplicitSurface() = default;
	virtual ~ImplicitSurface() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const { return std::abs(signedDistance(pos)); }
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool isInside(const VectorDr &pos) const { return isInside(signedDistance(pos)); }
};

template <int Dim>
class ImplicitSphere final : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDr _center;
	const real _radius;

public:

	ImplicitSphere(const VectorDr &center, const real radius) : _center(center), _radius(radius) { }
	virtual ~ImplicitSphere() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const override { return _center + closestNormal(pos) * _radius; }
	virtual VectorDr closestNormal(const VectorDr &pos) const override { return (pos - _center).normalized(); }
	virtual real signedDistance(const VectorDr &pos) const override { return (pos - _center).norm() - _radius; }
};

template <int Dim>
class ImplicitBox final : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDr _center;
	const VectorDr _halfLengths;

public:

	ImplicitBox(const VectorDr &minCorner, const VectorDr &lengths) : _center(minCorner + lengths / 2), _halfLengths(lengths / 2) { }
	virtual ~ImplicitBox() = default;

	virtual VectorDr closestNormal(const VectorDr &pos) const override
	{
		const VectorDr phi = (pos - _center).cwiseAbs() - _halfLengths;
		VectorDr normal;
		if ((phi.array() <= 0).all()) {
			int axis;
			phi.maxCoeff(&axis);
			normal = VectorDr::Unit(axis);
		}
		else normal = phi.cwiseMax(0);
		return normal.cwiseProduct((pos - _center).cwiseSign()).normalized();
	}

	virtual real signedDistance(const VectorDr &pos) const override
	{
		const VectorDr phi = (pos - _center).cwiseAbs() - _halfLengths;
		if ((phi.array() <= 0).all()) return phi.maxCoeff();
		else return phi.cwiseMax(0).norm();
	}
};

template <int Dim>
class ImplicitPlane final : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDr _position;
	const VectorDr _normal;

public:

	ImplicitPlane(const VectorDr &position, const VectorDr &direction) : _position(position), _normal(direction.normalized()) { }
	virtual ~ImplicitPlane() = default;

	virtual VectorDr closestNormal(const VectorDr &pos) const override { return _normal; }
	virtual real signedDistance(const VectorDr &pos) const override { return (pos - _position).dot(_normal); }
};

template <int Dim>
class ImplicitEllipsoid final : public ImplicitSurface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	const VectorDr _center;
	const VectorDr _semiAxels;

public:

	using Surface<Dim>::isInside;

	ImplicitEllipsoid(const VectorDr &center, const VectorDr &semiAxels) : _center(center), _semiAxels(semiAxels) { }
	virtual ~ImplicitEllipsoid() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const override { return pos / pos.cwiseQuotient(_semiAxels).norm(); } // not accurate solution
	virtual VectorDr closestNormal(const VectorDr &pos) const override { return (pos - closestPosition(pos)).normalized() * (isInside(pos) ? -1 : 1); }
	virtual real distance(const VectorDr &pos) const override { return (pos - closestPosition(pos)).norm(); }
	virtual real signedDistance(const VectorDr &pos) const override { return distance(pos) * (isInside(pos) ? -1 : 1); }
	virtual bool isInside(const VectorDr &pos) const override { return isInside(pos.cwiseQuotient(_semiAxels).squaredNorm() - 1); }
};

}
