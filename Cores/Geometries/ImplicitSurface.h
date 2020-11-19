#pragma once

#include "Geometries/Surface.h"

#include <cmath>

namespace PhysX {

template <int Dim>
class ImplicitSurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

public:

	using Surface<Dim>::inside;

	ImplicitSurface() = default;
	virtual ~ImplicitSurface() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const { return pos - signedDistance(pos) * closestNormal(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const { return std::abs(signedDistance(pos)); }
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool inside(const VectorDr &pos) const { return inside(signedDistance(pos)); }
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

	ImplicitBox(const VectorDr &center, const VectorDr &halfLengths) : _center(center), _halfLengths(halfLengths) { }
	virtual ~ImplicitBox() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const override
	{
		const VectorDr phi = (pos - _center).cwiseAbs() - _halfLengths;
		VectorDr sign = (pos - _center).cwiseSign();
		if ((phi.array() <= 0).all()) sign = sign.cwiseProduct(phi.cwiseEqual(phi.maxCoeff()));
		else sign = sign.cwiseProduct(phi.cwiseMax(0).cwiseSign());
		return _center + sign.cwiseProduct(_halfLengths);
	}

	virtual VectorDr closestNormal(const VectorDr &pos) const override
	{
		const VectorDr phi = (pos - _center).cwiseAbs() - _halfLengths;
		VectorDr sign = (pos - _center).cwiseSign();
		if ((phi.array() <= 0).all()) sign = sign.cwiseProduct(phi.cwiseEqual(phi.maxCoeff()));
		else sign = sign.cwiseProduct(phi.cwiseMax(0).cwiseSign());
		return sign.cwiseProduct(phi).normalized();
	}

	virtual real signedDistance(const VectorDr &pos) const override
	{
		const VectorDr phi = (pos - _center).cwiseAbs() - _halfLengths;
		if ((phi.array() <= 0).all()) return phi.maxCoeff();
		else return phi.cwiseMax(0).norm();
	}
};

}
