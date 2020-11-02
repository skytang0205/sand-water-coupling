#include "GeometryGenerator.h"

#include <iostream>

namespace PhysX {

GeometryGenerator::Data GeometryGenerator::createBox(float width, float height, float depth)
{
	Data ret;
	Vector3f halfLength = Vector3f(width, height, depth) * 0.5f;

	Vector3i coeff;
	for (int axis = 0; axis < 3; axis++) {
		for (int dir = 0; dir < 2; dir++) {
			uint base = uint(ret.vertices.size());
			coeff[axis] = dir * 2 - 1;
			for (int i = 0; i < 4; i++) {
				coeff[(axis + 1) % 3] = (i & 1) * 2 - 1;
				coeff[(axis + 2) % 3] = (i & 2) - 1;
				Vector3f pos = halfLength.cwiseProduct(coeff.cast<float>());
				Vector3f normal = Vector3f::Unit(axis) * coeff[axis];
				ret.vertices.push_back({ pos, normal });
			}
			if (dir) {
				ret.indices32.insert(ret.indices32.end(), { base + 0, base + 1, base + 3 });
				ret.indices32.insert(ret.indices32.end(), { base + 3, base + 2, base + 0 });
			}
			else {
				ret.indices32.insert(ret.indices32.end(), { base + 3, base + 1, base + 0 });
				ret.indices32.insert(ret.indices32.end(), { base + 0, base + 2, base + 3 });
			}
		}
	}
	return ret;
}

}
