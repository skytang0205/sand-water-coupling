#pragma once

#include "Graphics/GlRenderItem.h"
#include "Utilities/Yaml.h"

#include <string>
#include <unordered_map>
#include <vector>

namespace PhysX {

class GlSimulated : public GlRenderItem
{
public:

	static inline std::unordered_map<std::string, GLenum> strToMode = {
			{ "point_list", GL_POINTS },
			{ "line_list", GL_LINES },
			{ "triangle_list", GL_TRIANGLES }
		};

protected:

	Vector4f _diffuseAlbedo;
	Vector3f _fresnelR0;
	float _roughness;
	bool _enableColorMap;

	uint _frame = 0;

	std::vector<uint> _vtxFrameOffset;
	std::vector<uint> _idxFrameOffset;
	std::vector<Matrix4f> _worldMats;

public:

	GlSimulated(GlProgram *program, const std::string &outputDir, const uint endFrame, const int dim, const YAML::Node &node);

	GlSimulated() = delete;
	GlSimulated(const GlSimulated &rhs) = delete;
	GlSimulated &operator=(const GlSimulated &rhs) = delete;
	virtual ~GlSimulated() { }

	virtual void beginDraw() const override;
	void setFrame(const uint frame) { _frame = frame; }
	bool isTransparent() const { return _diffuseAlbedo.w() < 1.0f; }

protected:

	void loadMesh(const std::string &fileName, const int dim);
};

}
