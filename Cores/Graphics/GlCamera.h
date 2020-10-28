#pragma once

#include <Utilities/Types.h>

#include <algorithm>
#include <numbers>

#include <cmath>

namespace PhysX {

class GlCamera
{

protected:

	float _fovy;
	float _aspect;
	float _zNear;
	float _zFar;

	float _yScale = 1.0f;
	float _xScale = 1.0f;
	
	Vector3f _pos;
	Vector3f _target;
	Vector3f _worldUp;

	Vector3f _front = Vector3f::Zero();
	Vector3f _up = Vector3f::Zero();
	Vector3f _right = Vector3f::Zero();

	Matrix4f _proj = Matrix4f::Identity();
	Matrix4f _view = Matrix4f::Identity();
	Matrix4f _projView = Matrix4f::Identity();

	bool _projDirty;
	bool _viewDirty;

public:

	GlCamera(
		const float fovy,
		const float aspect,
		const float zNear,
		const float zFar,
		const Vector3f &pos,
		const Vector3f &target,
		const Vector3f &worldUp = Vector3f::Unit(1))
	{
		setPerspective(fovy, aspect, zNear, zFar);
		lookAt(pos, target, worldUp);
	}

	GlCamera() = delete;
	GlCamera(const GlCamera &rhs) = default;
	GlCamera &operator=(const GlCamera &rhs) = default;
	virtual ~GlCamera() = default;

	void update()
	{
		if (updateProjMatrix() || updateViewMatrix())
			_projView = _proj * _view;
	}

	void setPerspective(const float fovy, const float aspect, const float zNear, const float zFar)
	{
		_fovy = fovy;
		_aspect = aspect;
		_zNear = zNear;
		_zFar = zFar;
		_projDirty = true;
	}

	void setAspectRatio(const float aspect)
	{
		_aspect = aspect;
		_projDirty = true;
	}

	void lookAt(const Vector3f &pos, const Vector3f &target, const Vector3f &worldUp = Vector3f::Unit(1))
	{
		_pos = pos;
		_target = target;
		_worldUp = worldUp;
		_viewDirty = true;
	}

protected:

	virtual bool updateProjMatrix()
	{
		if (_projDirty) {
			_yScale = 1.0f / std::tan(_fovy / 2);
			_xScale = _yScale / _aspect;

			_proj << _xScale, 0, 0, 0,
				0, _yScale, 0, 0,
				0, 0, -(_zFar + _zNear) / (_zFar - _zNear), -1,
				0, 0, -2 * _zNear * _zFar / (_zFar - _zNear), 0;

			return _projDirty = false, true;
		}
		else return false;
	}

	virtual bool updateViewMatrix()
	{
		if (_viewDirty) {
			_front = (_target - _pos).normalized();
			_right = _front.cross(_worldUp).normalized();
			_up = _right.cross(_front).normalized();

			_view << _right, 0, _up, 0, -_front, 0, -_right.dot(_pos), -_up.dot(_pos), _front.dot(_pos);

			return _viewDirty = false, true;
		}
		else return false;
	}
};

class GlOrbitCamera : public GlCamera
{
protected:

	float _phi;
	float _theta;
	float _radius;

public:

	GlOrbitCamera(
		const float fovy,
		const float aspect,
		const float zNear,
		const float zFar,
		const float phi,
		const float theta,
		const float radius,
		const Vector3f &target,
		const Vector3f &worldUp = Vector3f::Unit(1)
		) :
		GlCamera(fovy, aspect, zNear, zFar, sphericalToCartesian(_phi, _theta, _radius), target, worldUp)
	{
	}

	void rotate(const float dx, const float dy)
	{
		_phi -= 0.25f * dx;
		_theta -= 0.25f * dy;
		std::clamp(_phi, 0.1f, float(std::numbers::pi) - 0.1f);
		_viewDirty = true;
	}

	void translate(const float dx, const float dy)
	{
	}

protected:

	virtual bool updateViewMatrix() override
	{
		if (_viewDirty) {
			_pos = sphericalToCartesian(_theta, _phi, _radius);
			return GlCamera::updateViewMatrix();
		}
		else return false;
	}

	static Vector3f sphericalToCartesian(const float phi, const float theta, const float radius)
	{
		return Vector3f(std::sin(theta) * std::sin(phi), std::cos(theta), std::sin(theta) * std::cos(phi)) * radius;
	}
};

}
