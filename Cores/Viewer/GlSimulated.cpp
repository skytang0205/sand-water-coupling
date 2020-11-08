#include "GlSimulated.h"

#include "Utilities/IO.h"

#include <fmt/core.h>

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
	_roughness = dim > 2 ? 0.75f : 1.03125;
	if (node["roughness"]) _roughness = node["roughness"].as<float>();

	// Initialize meshes.
	std::vector<Vector3f> positions;
	std::vector<Vector3f> normals;
	std::vector<float> heats;
	std::vector<uint> indices;
	if (dataMode == "static") {
		_vtxFrameOffset.push_back(0);
		if (_indexed) _idxFrameOffset.push_back(0);
		_frameWorlds.push_back(Matrix4f::Identity());
		loadMesh(fmt::format("{}/0/{}.mesh", outputDir, name), dim, positions, normals, heats, indices);
	}
	else if (dataMode == "dynamic") {
		_vtxFrameOffset.push_back(0);
		if (_indexed) _idxFrameOffset.push_back(0);
		_frameWorlds.push_back(Matrix4f::Identity());
		for (uint frame = 0; frame < endFrame; frame++) {
			loadMesh(fmt::format("{}/{}/{}.mesh", outputDir, frame, name), dim, positions, normals, heats, indices);
		}
	}
	for (size_t frame = 1; frame < _vtxFrameOffset.size(); frame++) {
		_vtxFrameOffset[frame] += _vtxFrameOffset[frame - 1];
		if (_indexed) _idxFrameOffset[frame] += _idxFrameOffset[frame - 1];
	}

	// Create vertex buffers.
	auto sizeBase = positions.size() * sizeof(float);
	glCreateBuffers(1, &_vbo);
	glNamedBufferData(_vbo, (_enableColorMap ? 7 : 6) * sizeBase, nullptr, GL_STATIC_DRAW);		
	// Attribute 0.
	glNamedBufferSubData(_vbo, 0, 3 * sizeBase, positions.data());
	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 3 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 0, 0);
	glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glEnableVertexArrayAttrib(_vao, 0);
	// Attribute 1.
	glNamedBufferSubData(_vbo, 3 * sizeBase, 3 * sizeBase, normals.data());
	glVertexArrayVertexBuffer(_vao, 1, _vbo, 3 * sizeBase, 3 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 1, 1);
	glVertexArrayAttribFormat(_vao, 1, 3, GL_FLOAT, GL_FALSE, 0);
	glEnableVertexArrayAttrib(_vao, 1);
	// Attribute 2.
	if (_enableColorMap) {
		glNamedBufferSubData(_vbo, 6 * sizeBase, sizeBase, heats.data());
		glVertexArrayVertexBuffer(_vao, 2, _vbo, 6 * sizeBase, sizeof(float));
		glVertexArrayAttribBinding(_vao, 2, 2);
		glVertexArrayAttribFormat(_vao, 2, 1, GL_FLOAT, GL_FALSE, 0);
		glEnableVertexArrayAttrib(_vao, 2);
	}

	// Create index buffers.
	if (_indexed) {
		glCreateBuffers(1, &_ebo);
		glNamedBufferStorage(_ebo, indices.size() * sizeof(indices[0]), indices.data(), 0);
		glVertexArrayElementBuffer(_vao, _ebo);
	}
}

void GlSimulated::beginDraw()
{
	// Set parameters according to current frame.
	if (_indexed) {
		_count = _idxFrameOffset[_currentFrame + size_t(1)] - _idxFrameOffset[_currentFrame];
		_indices = reinterpret_cast<const void *>(_idxFrameOffset[_currentFrame] * sizeof(uint));
		_baseVertex = _vtxFrameOffset[_currentFrame];
	}
	else {
		_first = _vtxFrameOffset[_currentFrame];
		_count = _vtxFrameOffset[_currentFrame + size_t(1)] - _vtxFrameOffset[_currentFrame];
	}
	// Set uniforms.
	_program->setUniform("uWorld", _frameWorlds[_currentFrame < _frameWorlds.size() ? _currentFrame : 0]);
	_program->setUniform("uDiffuseAlbedo", _diffuseAlbedo);
	_program->setUniform("uFresnelR0", _fresnelR0);
	_program->setUniform("uRoughness", _roughness);
	_program->setUniform("uEnableColorMap", uint(_enableColorMap));
}

void GlSimulated::loadMesh(
	const std::string &fileName,
	const int dim,
	std::vector<Vector3f> &positions,
	std::vector<Vector3f> &normals,
	std::vector<float> &heats,
	std::vector<uint> &indices
	)
{
	std::ifstream fin(fileName, std::ios::binary);
	uint vtxCnt;
	IO::readValue(fin, vtxCnt);
	_vtxFrameOffset.push_back(vtxCnt);
	// Read positions.
	positions.resize(positions.size() + vtxCnt, Vector3f::Zero().eval());
	if (dim > 2)
		IO::readArray(fin, positions.data() + positions.size() - vtxCnt, vtxCnt);
	else {
		for (uint i = 0; i < vtxCnt; i++) {
			IO::read(fin, positions.data() + positions.size() - vtxCnt + i, sizeof(Vector2f));
		}
	}
	// Read normals.
	normals.resize(normals.size() + vtxCnt, Vector3f::Zero().eval());
	if (dim > 2)
		IO::readArray(fin, normals.data() + normals.size() - vtxCnt, vtxCnt);
	else {
		for (uint i = 0; i < vtxCnt; i++)
			normals[normals.size() - vtxCnt + i].z() = 1.0f;
	}
	// Read heats.
	if (_enableColorMap) {
		heats.resize(heats.size() + vtxCnt);
		IO::readArray(fin, heats.data() + heats.size() - vtxCnt, vtxCnt);
	}
	// Read indices.
	if (_indexed) {
		uint idxCnt;
		IO::readValue(fin, idxCnt);
		_idxFrameOffset.push_back(idxCnt);
		indices.resize(indices.size() + idxCnt);
		IO::readArray(fin, indices.data() + indices.size() - idxCnt, idxCnt);
	}
}

}
