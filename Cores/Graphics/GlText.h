#pragma once

#include "Graphics/BitmapFont.h"
#include "Graphics/GlRenderItem.h"

namespace PhysX {

class GlText : public GlRenderItem
{
protected:

	static inline BitmapFont _bmConsolas = {
			64,
			51,
			512,
			512,
			1,
			128,
			reinterpret_cast<const char *>(BitmapFont::kBitmapConsolas[0]),
			BitmapFont::kBitmapConsolasChars
		};

public:

	GlText(GlProgram *program) : GlRenderItem(program) { }
};


}
