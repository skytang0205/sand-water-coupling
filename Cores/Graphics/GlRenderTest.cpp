#include "GlRenderTest.h"

#include "GeometryGenerator.h"

namespace PhysX {

GlRenderTest::GlRenderTest(GlProgram *program) : GlRenderItem(program)
{
	auto box = GeometryGenerator::createBox(1.0f, 2.0f, 3.0f);

	glCreateBuffers(1, &_vbo);
	glNamedBufferStorage(_vbo, box.vertices.size() * sizeof(box.vertices[0]), box.vertices.data(), 0);

	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 6 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribFormat(_vao, 1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float));
	glVertexArrayAttribBinding(_vao, 0, 0);
	glVertexArrayAttribBinding(_vao, 1, 0);
	glEnableVertexArrayAttrib(_vao, 0);
	glEnableVertexArrayAttrib(_vao, 1);

	glCreateBuffers(1, &_ebo);
	glNamedBufferStorage(_ebo, box.indices32.size() * sizeof(box.indices32[0]), box.indices32.data(), 0);
	glVertexArrayElementBuffer(_vao, _ebo);

	_indexed = true;
	_count = box.indices32.size();
}

void GlRenderTest::beginDraw() const
{
	_program->setUniform("uWorld", Matrix4f::Identity().eval());
	_program->setUniform("uDiffuseAlbedo", Vector4f(1, 0, 0, 1));
	_program->setUniform("uFresnelR0", Vector3f(0.125f, 0.125f, 0.125f));
	_program->setUniform("uRoughness", 0.02041f);
}

}
