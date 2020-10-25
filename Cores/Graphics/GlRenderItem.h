#pragma once

#include "GlProgram.h"

namespace PhysX {

class GlRenderItem
{
protected:

	GlProgram *_program;

public:

	GlRenderItem(GlProgram *program) : _program(program) { }

	GlRenderItem() = delete;
	GlRenderItem(const GlRenderItem &rhs) = delete;
	GlRenderItem &operator=(const GlRenderItem &rhs) = delete;
	virtual ~GlRenderItem() = default;

	virtual void draw() const = 0;
};

}
