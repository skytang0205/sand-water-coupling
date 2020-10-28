#pragma once

#include <Utilities/Types.h>

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

	const Vector3f _savedPos;
	const Vector3f _savedTarget;

	Vector3f _pos;
	Vector3f _target;

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
		const Vector3f &target
	) :
		_savedPos(pos),
		_savedTarget(target)
	{
		setPerspective(fovy, aspect, zNear, zFar);
		lookAt(pos, target);
	}

	GlCamera() = delete;
	GlCamera(const GlCamera &rhs) = default;
	GlCamera &operator=(const GlCamera &rhs) = default;
	virtual ~GlCamera() = default;

	void update()
	{
		if (_projDirty)
			updateProjMatrix();
		if (_viewDirty)
			updateViewMatrix();
		if (_projDirty || _viewDirty)
			_projView = _proj * _view;
		_projDirty = false;
		_viewDirty = false;
	}

	const Matrix4f &projView() const { return _projView; }

	virtual void reset() { lookAt(_savedPos, _savedTarget); }

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

	void lookAt(const Vector3f &pos, const Vector3f &target)
	{
		_pos = pos;
		_target = target;

		_front = (_target - _pos).normalized();
		_right = _front.cross(Vector3f::Unit(1)).normalized();
		_up = _right.cross(_front).normalized();

		_viewDirty = true;
	}

protected:

	virtual void updateProjMatrix()
	{
		_yScale = 1.0f / std::tan(_fovy / 2);
		_xScale = _yScale / _aspect;

		_proj << _xScale, 0, 0, 0,
			0, _yScale, 0, 0,
			0, 0, -(_zFar + _zNear) / (_zFar - _zNear), -1,
			0, 0, -2 * _zNear * _zFar / (_zFar - _zNear), 0;
	}

	virtual void updateViewMatrix()
	{
		_view << _right.x(), _up.x(), -_front.x(), -_right.dot(_pos),
			_right.y(), _up.y(), -_front.y(), -_up.dot(_pos),
			_right.z(), _up.z(), -_front.z(), _front.dot(_pos),
			0, 0, 0, 1;
	}
};
}
