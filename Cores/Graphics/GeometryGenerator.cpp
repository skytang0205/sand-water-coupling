#include "GeometryGenerator.h"

#include <iostream>
#include <numbers>

#include <cmath>
#include <cstdlib>

namespace PhysX {

GeometryGenerator::Data GeometryGenerator::createBox(const Vector3f &lengths)
{
	Data ret;
	Vector3f halfLengths = lengths * 0.5f;

	Vector3i coeff;
	for (int axis = 0; axis < 3; axis++) {
		for (int dir = 0; dir < 2; dir++) {
			uint base = uint(ret.vertices.size());
			coeff[axis] = dir * 2 - 1;
			for (int i = 0; i < 4; i++) {
				coeff[(axis + 1) % 3] = (i & 1) * 2 - 1;
				coeff[(axis + 2) % 3] = (i & 2) - 1;
				Vector3f pos = halfLengths.cwiseProduct(coeff.cast<float>());
				Vector3f normal = Vector3f::Unit(axis) * coeff[axis];
				ret.vertices.push_back({ pos, 0, normal, 0 });
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

GeometryGenerator::Data GeometryGenerator::createUVSphere(const float radius, const uint sliceCnt, const uint stackCnt)
{
	if (radius <= 0.0f || sliceCnt < 3 || stackCnt < 2) {
		std::cerr << "Error: [GeometryGenerator] created UV sphere width invalid parameters." << std::endl;
		std::exit(-1);
	}
	Data ret;
	const float deltaPhi = float(std::numbers::pi) / stackCnt; // latitude, in [-pi / 2, +pi / 2]
	const float deltaLambda = float(std::numbers::pi) * 2.0f / sliceCnt; // longitude, in [-pi, pi]

	ret.vertices.push_back({ Vector3f::Unit(1) * (-radius), 0, Vector3f::Unit(1) * (-1), 0 });
	for (uint i = 1; i <= stackCnt - 1; i++) {
		const float phi = i * deltaPhi - float(std::numbers::pi) * 0.5f;
		for (uint j = 0; j < sliceCnt; j++) {
			const float lambda = j * deltaLambda - float(std::numbers::pi);
			Vertex vert;
			vert.normal = Vector3f(std::cos(phi) * std::sin(lambda), std::sin(phi), std::cos(phi) * std::cos(lambda));
			vert.enableColorMap = 0;
			vert.pos = vert.normal * radius;
			vert.heat = 0;
			ret.vertices.push_back(vert);
		}
	}
	ret.vertices.push_back({ Vector3f::Unit(1) * radius, 0, Vector3f::Unit(1), 0 });

	for (uint i = 1; i <= sliceCnt; i++)
		ret.indices32.insert(ret.indices32.end(), { 0, 1 + i % sliceCnt, i });
	for (uint i = 0; i < stackCnt - 2; i++)
		for (uint j = 0; j < sliceCnt; j++) {
			const uint index0 = 1 + i * sliceCnt + j;
			const uint index1 = 1 + i * sliceCnt + (j + 1) % sliceCnt;
			const uint index2 = 1 + (i + 1) * sliceCnt + j;
			const uint index3 = 1 + (i + 1) * sliceCnt + (j + 1) % sliceCnt;
			ret.indices32.insert(ret.indices32.end(), { index2, index0, index1 });
			ret.indices32.insert(ret.indices32.end(), { index1, index3, index2 });
		}
	const uint topIndex = uint(ret.vertices.size()) - 1;
	for (uint i = 1; i <= sliceCnt; i++)
		ret.indices32.insert(ret.indices32.end(), { topIndex, topIndex - i % sliceCnt - 1, topIndex - i });

	return ret;
}

}
