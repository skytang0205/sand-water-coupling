#include "GlRenderTest.h"

#include "GeometryGenerator.h"

namespace PhysX {

GlRenderTest::GlRenderTest(GlProgram *program) : GlRenderItem(program)
{
	auto box = GeometryGenerator::createBox(1.0f, 2.0f, 3.0f);
	static constexpr float colors[] = {
		1.0f, 0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 0.0f, 1.0f,
		0.0f, 1.0f, 0.0f, 1.0f,
		0.0f, 0.0f, 1.0f, 1.0f,
		0.0f, 0.0f, 1.0f, 1.0f,
		0.0f, 0.0f, 1.0f, 1.0f,
		0.0f, 0.0f, 1.0f, 1.0f,
		1.0f, 1.0f, 0.0f, 1.0f,
		1.0f, 1.0f, 0.0f, 1.0f,
		1.0f, 1.0f, 0.0f, 1.0f,
		1.0f, 1.0f, 0.0f, 1.0f,
		1.0f, 0.0f, 1.0f, 1.0f,
		1.0f, 0.0f, 1.0f, 1.0f,
		1.0f, 0.0f, 1.0f, 1.0f,
		1.0f, 0.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 1.0f, 1.0f,
		0.0f, 1.0f, 1.0f, 1.0f
	};
	std::vector<float> vertices;
	for (int i = 0; i < 24; i++) {
		vertices.insert(vertices.end(), box.vertices[i].pos.data(), box.vertices[i].pos.data() + 3);
		vertices.insert(vertices.end(), colors + i * 4, colors + i * 4 + 4);
	}

	glCreateBuffers(1, &_vbo);
	glNamedBufferStorage(_vbo, vertices.size() * sizeof(vertices[0]), vertices.data(), 0);

	glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 7 * sizeof(float));
	glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
	glVertexArrayAttribFormat(_vao, 1, 4, GL_FLOAT, GL_FALSE, 3 * sizeof(float));
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

}
