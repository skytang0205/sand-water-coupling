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

	GlRenderItem(GlProgram *program) : _program(program) { glGenVertexArrays(1, &_vao); }

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
		glGenBuffers(1, &_vbo);
		glGenBuffers(1, &_ebo);

		glBindVertexArray(_vao);

		const float vertices[] = {
			 0.5f,  0.5f, 0.0f,
			 0.5f, -0.5f, 0.0f,
			-0.5f, -0.5f, 0.0f,
			-0.5f,  0.5f, 0.0f 
		};
		const uint indices[] = {
			0, 1, 3,
			1, 2, 3
		};

		glBindBuffer(GL_ARRAY_BUFFER, _vbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _ebo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void *)0);
		glEnableVertexAttribArray(0);

		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		_indexed = true;
		_count = 6 / 2;
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
		_program->setUniform("uColor", Vector4f(1, 0, 0, 0));
	}

};

}
