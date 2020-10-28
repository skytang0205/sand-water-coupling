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

	float _yScale;
	float _xScale;
	
	Vector3f _pos;
	Vector3f _target;
	Vector3f _worldUp = Vector3f::Unit(1);

	Vector3f _front;
	Vector3f _up;
	Vector3f _right;

	Matrix4f _proj = Matrix4f::Identity();
	Matrix4f _view = Matrix4f::Identity();
	Matrix4f _projView;

	bool _projDirty = true;
	bool _viewDirty = true;

public:

	void update()
	{
		if (_projDirty) updateProjMatrix();
		if (_viewDirty) updateViewMatrix();
		if (_projDirty || _viewDirty) _projView = _proj * _view;
		_projDirty = false;
		_viewDirty = false;
	}

protected:

	void updateProjMatrix()
	{
		_yScale = 1.0f / std::tan(_fovy / 2);
		_xScale = _yScale / _aspect;

		_proj << _xScale, 0, 0, 0,
			0, _yScale, 0, 0,
			0, 0, -(_zFar + _zNear) / (_zFar - _zNear), -1,
			0, 0, -2 * _zNear * _zFar / (_zFar - _zNear), 0;
	}

	void updateViewMatrix()
	{
		_front = (_target - _pos).normalized();
		_right = _front.cross(_worldUp).normalized();
		_up = _right.cross(_front).normalized();

		_view << _right, 0, _up, 0, -_front, 0, -_right.dot(_pos), -_up.dot(_pos), _front.dot(_pos);
	}
};

}
