#pragma once

#include "Graphics/GlRenderItem.h"

namespace PhysX {

class GlRenderTest: public GlRenderItem
{
protected:

	GLuint _vbo;
	GLuint _ebo;

public:

	GlRenderTest(GlProgram *program);

	GlRenderTest(const GlRenderTest &rhs) = delete;
	GlRenderTest &operator=(const GlRenderTest &rhs) = delete;

	virtual ~GlRenderTest()
	{
		glDeleteBuffers(1, &_ebo);
		glDeleteBuffers(1, &_vbo);
	}

	virtual void beginDraw() override;
};

}
