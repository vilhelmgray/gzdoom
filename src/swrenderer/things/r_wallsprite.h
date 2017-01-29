//
// Copyright (C) 1993-1996 by id Software, Inc.
//
// This source is available for distribution and/or modification
// only under the terms of the DOOM Source Code License as
// published by id Software. All rights reserved.
//
// The source is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// FITNESS FOR A PARTICULAR PURPOSE. See the DOOM Source Code License
// for more details.
//

#pragma once

#include "r_visiblesprite.h"

namespace swrenderer
{
	class ProjectedWallTexcoords;
	class ColumnDrawerArgs;

	class RenderWallSprite : public VisibleSprite
	{
	public:
		static void Project(AActor *thing, const DVector3 &pos, FTextureID picnum, const DVector2 &scale, int renderflags, int spriteshade, bool foggy, FDynamicColormap *basecolormap);

	protected:
		bool IsWallSprite() const override { return true; }
		void Render(short *cliptop, short *clipbottom, int minZ, int maxZ) override;

	private:
		static void DrawColumn(ColumnDrawerArgs &drawerargs, int x, FTexture *WallSpriteTile, const ProjectedWallTexcoords &walltexcoords, double texturemid, float maskedScaleY, bool sprflipvert, const short *mfloorclip, const short *mceilingclip);

		FWallCoords wallc;
		uint32_t Translation = 0;
		uint32_t FillColor = 0;
	};
}
