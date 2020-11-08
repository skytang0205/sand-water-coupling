#include "GlSimulated.h"

#include <fmt/core.h>

#include <iostream>

#include <cstdlib>

namespace PhysX {
GlSimulated::GlSimulated(GlProgram *program, const std::string &outputDir, const uint endFrame, const int dim, const YAML::Node &node) :
	GlRenderItem(program)
{
	// Read basic properties.
	const std::string name = node["name"].as<std::string>();
	const std::string dataMode = node["data_mode"].as<std::string>();
	const std::string primitiveType = node["primitive_type"].as<std::string>();
	if (name.empty() || dataMode.empty() || strToMode.find(primitiveType) == strToMode.end()) {
		std::cerr << fmt::format("Error: [GlSimulated] encountered invalid simulated object \"{}\".", name) << std::endl;
		std::exit(-1);
	}
	_mode = strToMode[primitiveType];
	_indexed = node["indexed"].as<bool>();

	// Read properties used in shader.
	_enableColorMap = node["enable_color_map"].as<bool>();
	_diffuseAlbedo = Vector4f(0.5f, 0.5f, 0.5f, 1.0f);
	if (node["diffuse_albedo"]) _diffuseAlbedo = node["diffuse_albedo"].as<Vector4f>();
	_fresnelR0 = Vector3f(0.02041f, 0.02041f, 0.02041f);
	if (node["fresnel_r0"]) _fresnelR0 = node["fresnel_r0"].as<Vector3f>();
	_roughness = dim > 2 ? 0.75f : 1.125f;
	if (node["roughness"]) _roughness = node["roughness"].as<float>();
}

void GlSimulated::beginDraw() const
{
	_program->setUniform("uDiffuseAlbedo", _diffuseAlbedo);
	_program->setUniform("uFresnelR0", _fresnelR0);
	_program->setUniform("uRoughness", _roughness);
	_program->setUniform("uEnableColorMap", uint(_enableColorMap));
}

}
