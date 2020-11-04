#include "GlRenderTest.h"

#include "GeometryGenerator.h"

namespace PhysX {

GlRenderTest::GlRenderTest(GlProgram *program) : GlRenderItem(program)
{
	auto geo = GeometryGenerator::createBox(1.0, 1.0, 1.0);

	glCreateBuffers(1, &_vbo);
	glNamedBufferStorage(_vbo, geo.vertices.size() * sizeof(geo.vertices[0]), geo.vertices.data(), 0);

	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 8 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribFormat(_vao, 1, 1, GL_UNSIGNED_INT, GL_FALSE, 3 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 2, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 3, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 0, 0);
	glVertexArrayAttribBinding(_vao, 1, 0);
	glVertexArrayAttribBinding(_vao, 2, 0);
	glVertexArrayAttribBinding(_vao, 3, 0);
	glEnableVertexArrayAttrib(_vao, 0);
	glEnableVertexArrayAttrib(_vao, 1);
	glEnableVertexArrayAttrib(_vao, 2);
	glEnableVertexArrayAttrib(_vao, 3);

	glCreateBuffers(1, &_ebo);
	glNamedBufferStorage(_ebo, geo.indices32.size() * sizeof(geo.indices32[0]), geo.indices32.data(), 0);
	glVertexArrayElementBuffer(_vao, _ebo);

	_indexed = true;
	_count = uint(geo.indices32.size());
}

void GlRenderTest::beginDraw() const
{
	_program->setUniform("uWorld", Matrix4f::Identity().eval());
	_program->setUniform("uDiffuseAlbedo", Vector4f(1, 0, 0, 1));
	_program->setUniform("uFresnelR0", (0.02041f * Vector3f::Ones()).eval());
	_program->setUniform("uRoughness", 0.75f);
}

}
