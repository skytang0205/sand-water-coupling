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
		Vector3f normal;
	};

	struct Data
	{
		std::vector<Vertex> vertices;
		std::vector<uint> indices32;
	};

public:

	static Data createBox(float width, float height, float depth);
};

}
