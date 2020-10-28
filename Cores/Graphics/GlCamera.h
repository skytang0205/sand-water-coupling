#pragma once

#include <Utilities/Types.h>

namespace PhysX {

class GlCamera
{

protected:

	Vector3f _pos;
	Vector3f _target;

	Vector3f _up;
	Vector3f _right;
	Vector3f _worldUp;

	Matrix4f _view;

public:

	virtual void update()
	{
	}

protected:

	static Matrix4f lookAt(const Vector3f &pos, const Vector3f &target, const Vector3f &up)
	{
		Vector3f front = target - pos;
	}
};

}
