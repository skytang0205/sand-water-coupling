#pragma once

namespace PhysX {

struct BitmapFont final
{
	struct CharDesc
	{
		int x;
		int y;
		int width;
		int height;
		int xOffset;
		int yOffset;
		int xAdvance;
	};

#include "BitmapConsolas.inc"

	int lineHeight;
	int base;
	int scaleW;
	int scaleH;
	int nChnls;
	int count;
	const char *data;
	const CharDesc *chars;
};

};
