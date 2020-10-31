#pragma once

namespace PhysX {

struct BitmapFont final
{
	struct CharDesc
	{
		int id;
		int x;
		int y;
		int width;
		int height;
		int xOffset;
		int yOffset;
		int xadvance;
	};

#include "BmConsolas.inc"

	int lineHeight;
	int base;
	int scaleW;
	int scaleH;
	int nChnls;
	const char *data;
	CharDesc *chars;
};

}
