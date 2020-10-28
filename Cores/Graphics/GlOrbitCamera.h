#pragma once

#include "Graphics/GlCamera.h"

#include <algorithm>
#include <numbers>

namespace PhysX {

class GlOrbitCamera : public GlCamera
{
protected:

	const float _savedRadius;
	const float _savedPhi;
	const float _savedTheta;

	float _radius;
	float _phi;
	float _theta;

public:

	GlOrbitCamera(
		const float fovy,
		const float aspect,
		const float zNear,
		const float zFar,
		const float radius,
		const float phi,
		const float theta,
		const Vector3f &target
	) :
		GlCamera(fovy, aspect, zNear, zFar, sphericalToCartesian(radius, phi, theta), target),
		_savedRadius(radius),
		_savedPhi(phi),
		_savedTheta(theta),
		_radius(radius),
		_phi(phi),
		_theta(theta)
	{ }

	GlOrbitCamera() = delete;
	GlOrbitCamera(const GlOrbitCamera &rhs) = default;
	GlOrbitCamera &operator=(const GlOrbitCamera &rhs) = default;
	virtual ~GlOrbitCamera() = default;

	virtual void reset() override
	{
		GlCamera::reset();
		setSpherical(_savedRadius, _savedPhi, _savedTheta);
	}

	void setSpherical(const float radius, const float phi, const float theta)
	{
		_radius = radius;
		_phi = phi;
		_theta = theta;
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	void rotate(const float dx, const float dy)
	{
		static constexpr float kRotateRatio = 0.25f * float(std::numbers::pi) / 180.0f;
		_phi += kRotateRatio * dx;
		_theta += kRotateRatio * dy;
		_theta = std::clamp(_theta, 0.1f, float(std::numbers::pi) - 0.1f);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	void translate(const float dx, const float dy)
	{
		static constexpr float kTranslateRatio = 0.001f;
		_target += kTranslateRatio * _radius * (dy * _up - dx * _right);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

	void scale(const float dy)
	{
		static constexpr float kScaleRatio = 0.1f;
		_radius /= std::exp(kScaleRatio * dy);
		_radius = std::clamp(_radius, 0.1f, 150.0f);
		lookAt(sphericalToCartesian(_radius, _phi, _theta), _target);
	}

protected:

	static Vector3f sphericalToCartesian(const float radius, const float phi, const float theta)
	{
		return Vector3f(std::sin(theta) * std::sin(phi), std::cos(theta), std::sin(theta) * std::cos(phi)) * radius;
	}
};

}
