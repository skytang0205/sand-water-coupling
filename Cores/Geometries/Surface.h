#pragma once

#include "Utilities/Types.h"

#include <memory>

namespace PhysX {

template <int Dim>
class Surface
{
	DECLARE_DIM_TYPES(Dim)

public:

	Surface() = default;
	virtual ~Surface() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const = 0;
	virtual VectorDr closestNormal(const VectorDr &pos) const = 0;
	virtual real distance(const VectorDr &pos) const = 0;
	virtual real signedDistance(const VectorDr &pos) const = 0;
	virtual bool isInside(const VectorDr &pos) const = 0;

	static bool isInside(const real phi) { return phi <= 0; }
	static int sign(const real phi) { return isInside(phi) ? -1 : 1; }
	static bool isInterface(const real phi0, const real phi1) { return sign(phi0) != sign(phi1); }
	static real theta(const real phi0, const real phi1) { return phi0 / (phi0 - phi1); }
	static real fraction(const real phi0, const real phi1)
	{
		if (isInside(phi0) && isInside(phi1)) return 1;
		else if (isInside(phi0) && !isInside(phi1)) return theta(phi0, phi1);
		else if (!isInside(phi0) && isInside(phi1)) return theta(phi1, phi0);
		else return 0;
	}
};

template <int Dim>
class ComplementarySurface : public Surface<Dim>
{
	DECLARE_DIM_TYPES(Dim)

protected:

	std::unique_ptr<Surface<Dim>> _surface;

public:

	ComplementarySurface(std::unique_ptr<Surface<Dim>> surface) : _surface(std::move(surface)) { }
	virtual ~ComplementarySurface() = default;

	virtual VectorDr closestPosition(const VectorDr &pos) const override { return _surface->closestPosition(pos); }
	virtual VectorDr closestNormal(const VectorDr &pos) const override { return -_surface->closestNormal(pos); }
	virtual real distance(const VectorDr &pos) const override { return _surface->distance(pos); }
	virtual real signedDistance(const VectorDr &pos) const override { return -_surface->signedDistance(pos); }
	virtual bool isInside(const VectorDr &pos) const override { return Surface<Dim>::isInside(signedDistance(pos)); }
};

}
