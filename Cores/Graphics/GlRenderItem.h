#pragma once

#include "Graphics/GlProgram.h"

namespace PhysX {

class GlRenderItem
{
protected:

	GLenum _mode = GL_TRIANGLES;
	GLint _first = 0;
	GLsizei _count = 0;
	GLenum _type = GL_UNSIGNED_INT;
	const void *_indices = 0;
	GLsizei _instanceCount = 1;
	GLint _baseVertex = 0;
	GLuint _baseInstance = 0;

	bool _indexed = false;
	bool _visible = true;

	GlProgram *_program;

	GLuint _vao;

public:

	GlRenderItem(GlProgram *program) : _program(program) { glCreateVertexArrays(1, &_vao); }

	GlRenderItem() = delete;
	GlRenderItem(const GlRenderItem &rhs) = delete;
	GlRenderItem &operator=(const GlRenderItem &rhs) = delete;
	virtual ~GlRenderItem() { glDeleteVertexArrays(1, &_vao); }

	virtual void beginDraw() const { }
	virtual void endDraw() const { }

	void draw() const;

	bool isVisible() const { return _visible; }
	void show() { _visible = true; }
	void hide() { _visible = false; }
};

class GlRenderTest : public GlRenderItem
{
protected:

	GLuint _vbo;
	GLuint _ebo;

public:

	GlRenderTest(GlProgram *program) : GlRenderItem(program)
	{
		static constexpr float vertices[] = {
			-0.5f, -0.5f, -0.5f, 0.0f, 0.0f, 0.0f, 0.5f,
			-0.5f, -0.5f,  0.5f, 0.0f, 0.0f, 1.0f, 0.5f,
			-0.5f,  0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.5f,
			-0.5f,  0.5f,  0.5f, 0.0f, 1.0f, 1.0f, 0.5f,
			 0.5f, -0.5f, -0.5f, 1.0f, 0.0f, 0.0f, 0.5f,
			 0.5f, -0.5f,  0.5f, 1.0f, 0.0f, 1.0f, 0.5f,
			 0.5f,  0.5f, -0.5f, 1.0f, 1.0f, 0.0f, 0.5f,
			 0.5f,  0.5f,  0.5f, 1.0f, 1.0f, 1.0f, 0.5f
		};
		static constexpr uint indices[] = {
			3, 2, 0,
			0, 1, 3,
			4, 6, 7,
			7, 5, 4,
			7, 3, 1,
			1, 5, 7,
			0, 2, 6,
			6, 4, 0,
			5, 1, 0,
			0, 4, 5,
			2, 3, 7,
			7, 6, 2
		};

		glCreateBuffers(1, &_vbo);
		glNamedBufferStorage(_vbo, sizeof(vertices), vertices, 0);

		glVertexArrayVertexBuffer(_vao, 0, _vbo, 0, 7 * sizeof(float));
		glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
		glVertexArrayAttribFormat(_vao, 1, 4, GL_FLOAT, GL_FALSE, 3 * sizeof(float));
		glVertexArrayAttribBinding(_vao, 0, 0);
		glVertexArrayAttribBinding(_vao, 1, 0);
		glEnableVertexArrayAttrib(_vao, 0);
		glEnableVertexArrayAttrib(_vao, 1);

		glCreateBuffers(1, &_ebo);
		glNamedBufferStorage(_ebo, sizeof(indices), indices, 0);
		glVertexArrayElementBuffer(_vao, _ebo);

		_indexed = true;
		_count = 36;
	}

	GlRenderTest() = delete;
	GlRenderTest(const GlRenderTest &rhs) = delete;
	GlRenderTest &operator=(const GlRenderTest &rhs) = delete;
	virtual ~GlRenderTest()
	{
		glDeleteBuffers(1, &_ebo);
		glDeleteBuffers(1, &_vbo);
	}

	virtual void beginDraw() const override
	{
	}

};

}
