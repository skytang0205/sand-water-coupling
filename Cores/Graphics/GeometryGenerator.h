#pragma once

#include "Utilities/Types.h"

#include <vector>

namespace PhysX {

class GeometryGenerator
{
public:

	struct Vertex
	{
		Vector3f pos;
		uint enableColorMap;
		Vector3f normal;
		float heat;
	};

	struct Data
	{
		std::vector<Vertex> vertices;
		std::vector<uint> indices32;
	};

public:

	static Data createBox(const Vector3f &lengths);
	static Data createBox(const float width, const float height, const float depth) { return createBox(Vector3f(width, height, depth)); }
	static Data createUVSphere(const float radius, const uint sliceCnt, const uint stackCnt);
};

}
