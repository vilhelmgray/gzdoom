/*
**  Polygon Doom software renderer
**  Copyright (c) 2016 Magnus Norddahl
**
**  This software is provided 'as-is', without any express or implied
**  warranty.  In no event will the authors be held liable for any damages
**  arising from the use of this software.
**
**  Permission is granted to anyone to use this software for any purpose,
**  including commercial applications, and to alter it and redistribute it
**  freely, subject to the following restrictions:
**
**  1. The origin of this software must not be misrepresented; you must not
**     claim that you wrote the original software. If you use this software
**     in a product, an acknowledgment in the product documentation would be
**     appreciated but is not required.
**  2. Altered source versions must be plainly marked as such, and must not be
**     misrepresented as being the original software.
**  3. This notice may not be removed or altered from any source distribution.
**
*/

#pragma once

#include "screen_triangle.h"
#include "swfragmentshader32.h"

class Sampler8
{
public:
	uint32_t width;
	uint32_t height;
	const uint8_t *data;

	FORCEINLINE uint8_t TextureNearest(const SWVec4f &input, int i) const;
};

class TranslatedSampler8
{
public:
	uint32_t width;
	uint32_t height;
	const uint8_t *data;
	const uint8_t *translation;

	FORCEINLINE uint8_t TextureNearest(const SWVec4f &input, int i) const;
};

template<typename ModeT>
class SWFragmentShader8
{
public:
	// Uniforms
	Sampler8 Tex;
	TranslatedSampler8 TranslatedTex;
	const uint8_t *Colormaps;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;
	PolyLight *Lights;
	int NumLights;
	float WorldNormal[3];
	uint8_t DynLightColor[4];
	uint32_t SkycapColor;
	uint8_t FillColor;
	uint32_t Alpha;

	// In variables
	float GradW[4];
	SWVec4f WorldPos;
	SWVec4f TexCoord;

	// Out variables
	uint8_t FragColor[4];
	uint8_t FragAlpha[4];

	FORCEINLINE void SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void Run();

private:
	FORCEINLINE uint32_t CalcDynamicLight(int i);
};

template<typename ModeT>
class ScreenBlockDrawer8
{
public:
	static void Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args);

private:
	FORCEINLINE void SetUniforms(const TriDrawTriangleArgs *args);
	FORCEINLINE void SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY);
	FORCEINLINE void ProcessBlock(uint32_t mask0, uint32_t mask1);
	FORCEINLINE void ProcessMaskRange(uint32_t mask);
	FORCEINLINE void StepY();
	FORCEINLINE void StoreFull(int offset);
	FORCEINLINE void StoreMasked(int offset, uint32_t mask);
	FORCEINLINE uint8_t Blend(uint8_t *destptr, uint8_t src, uint8_t srcalpha);

	// Gradients
	ScreenTriangleStepVariables GradPosX;
	ScreenTriangleStepVariables GradPosY;
	ScreenTriangleStepVariables GradStepX;
	ScreenTriangleStepVariables GradStepY;

	// Blend stage
	uint8_t *Dest;
	int Pitch;

	SWFragmentShader8<ModeT> Shader;
};

/////////////////////////////////////////////////////////////////////////////

uint8_t Sampler8::TextureNearest(const SWVec4f &input, int i) const
{
	uint32_t x = ((static_cast<uint32_t>(static_cast<int32_t>(input.x[i] * (1 << 24)) << 8) >> 16) * width) >> 16;
	uint32_t y = ((static_cast<uint32_t>(static_cast<int32_t>(input.y[i] * (1 << 24)) << 8) >> 16) * height) >> 16;
	return data[y + x * height];
}

/////////////////////////////////////////////////////////////////////////////

uint8_t TranslatedSampler8::TextureNearest(const SWVec4f &input, int i) const
{
	uint32_t x = ((static_cast<uint32_t>(static_cast<int32_t>(input.x[i] * (1 << 24)) << 8) >> 16) * width) >> 16;
	uint32_t y = ((static_cast<uint32_t>(static_cast<int32_t>(input.y[i] * (1 << 24)) << 8) >> 16) * height) >> 16;
	return translation[data[y + x * height]];
}

/////////////////////////////////////////////////////////////////////////////

template<typename ModeT>
void SWFragmentShader8<ModeT>::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
{
	for (int i = 0; i < 4; i++)
	{
		GradW[i] = pos.W;
		float w = 1.0f / pos.W;

		WorldPos.x[i] = pos.WorldX * w;
		WorldPos.y[i] = pos.WorldY * w;
		WorldPos.z[i] = pos.WorldZ * w;
		TexCoord.x[i] = pos.U * w;
		TexCoord.y[i] = pos.V * w;

		pos.W += step.W;

		pos.WorldX += step.WorldX;
		pos.WorldY += step.WorldY;
		pos.WorldZ += step.WorldZ;
		pos.U += step.U;
		pos.V += step.V;
	}
}

template<typename ModeT>
void SWFragmentShader8<ModeT>::Run()
{
	using namespace TriScreenDrawerModes;

	for (int i = 0; i < 4; i++)
	{
		uint8_t fg;
		uint8_t fgalpha = 255;

		if (ModeT::SWFlags & SWSTYLEF_Fill)
		{
			fg = FillColor;
		}
		else if (ModeT::SWFlags & SWSTYLEF_FogBoundary)
		{
			fg = FragColor[i];
		}
		else if (ModeT::SWFlags & SWSTYLEF_Translated)
		{
			fg = TranslatedTex.TextureNearest(TexCoord, i);
			fgalpha = (fg != 0) ? 255 : 0;
		}
		else
		{
			fg = Tex.TextureNearest(TexCoord, i);
			fgalpha = (fg != 0) ? 255 : 0;
		}

		if ((ModeT::Flags & STYLEF_ColorIsFixed) && !(ModeT::SWFlags & SWSTYLEF_Fill))
		{
			if (ModeT::Flags & STYLEF_RedIsAlpha)
				fgalpha = fg;
			fg = FillColor;
		}

		if (!(ModeT::Flags & STYLEF_Alpha1))
		{
			fgalpha = (fgalpha * Alpha) >> 8;
		}

		if (ModeT::SWFlags & SWSTYLEF_Skycap)
		{
			int32_t v = static_cast<int32_t>(TexCoord.y[i] * (1 << 24));

			int start_fade = 2; // How fast it should fade out

			int alpha_top = clamp(v >> (16 - start_fade), 0, 256);
			int alpha_bottom = clamp(((2 << 24) - v) >> (16 - start_fade), 0, 256);
			int a = MIN(alpha_top, alpha_bottom);
			int inv_a = 256 - a;

			if (a == 256)
			{
				FragColor[i] = fg;
			}
			else
			{
				uint32_t capcolor = GPalette.BaseColors[SkycapColor].d;
				uint32_t texelrgb = GPalette.BaseColors[fg].d;

				uint32_t r = RPART(texelrgb);
				uint32_t g = GPART(texelrgb);
				uint32_t b = BPART(texelrgb);
				uint32_t fg_a = APART(texelrgb);
				uint32_t bg_red = RPART(capcolor);
				uint32_t bg_green = GPART(capcolor);
				uint32_t bg_blue = BPART(capcolor);
				r = (r * a + bg_red * inv_a + 127) >> 8;
				g = (g * a + bg_green * inv_a + 127) >> 8;
				b = (b * a + bg_blue * inv_a + 127) >> 8;

				FragColor[i] = RGB256k.All[((r >> 2) << 12) | ((g >> 2) << 6) | (b >> 2)];
			}

			FragAlpha[i] = fgalpha;
		}
		else
		{
			float lightposf;
			if (!LightMask)
				lightposf = clamp(Shade - MIN(24.0f / 32.0f, GlobVis * GradW[i]), 0.0f, 31.0f / 32.0f);
			else
				lightposf = 1.0f - Light;

			uint32_t lightshade = ((int32_t)(lightposf * NUMCOLORMAPS)) << 8;
			uint8_t shadedfg = Colormaps[lightshade + fg];

			if (NumLights != 0)
			{
				uint32_t texelrgb = GPalette.BaseColors[fg].d;
				uint32_t material_r = RPART(texelrgb);
				uint32_t material_g = GPART(texelrgb);
				uint32_t material_b = BPART(texelrgb);

				uint32_t shadedrgb = GPalette.BaseColors[shadedfg].d;
				uint32_t r = RPART(shadedrgb);
				uint32_t g = GPART(shadedrgb);
				uint32_t b = BPART(shadedrgb);

				uint32_t dynlight = CalcDynamicLight(i);
				r = MIN(r + ((material_r * RPART(dynlight)) >> 8), (uint32_t)255);
				g = MIN(g + ((material_g * GPART(dynlight)) >> 8), (uint32_t)255);
				b = MIN(b + ((material_b * BPART(dynlight)) >> 8), (uint32_t)255);

				shadedfg = RGB256k.All[((r >> 2) << 12) | ((g >> 2) << 6) | (b >> 2)];
			}

			FragColor[i] = shadedfg;
			FragAlpha[i] = fgalpha;
		}
	}
}

template<typename ModeT>
uint32_t SWFragmentShader8<ModeT>::CalcDynamicLight(int j)
{
	uint32_t lit[4];

	lit[0] = DynLightColor[0];
	lit[1] = DynLightColor[1];
	lit[2] = DynLightColor[2];
	lit[3] = DynLightColor[3];

	for (int i = 0; i < NumLights; i++)
	{
		float light_radius = Lights[i].radius;

		bool is_attenuated = light_radius < 0.0f;
		if (is_attenuated)
			light_radius = -light_radius;

		// L = light-pos
		// dist = sqrt(dot(L, L))
		// distance_attenuation = 1 - MIN(dist * (1/radius), 1)
		float Lx = Lights[i].x - WorldPos.x[j];
		float Ly = Lights[i].y - WorldPos.y[j];
		float Lz = Lights[i].z - WorldPos.z[j];
		float dist2 = Lx * Lx + Ly * Ly + Lz * Lz;
		float rcp_dist = 1.0f / sqrt(dist2);
		float dist = dist2 * rcp_dist;
		float distance_attenuation = 256.0f - MIN(dist * light_radius, 256.0f);

		// The simple light type
		float simple_attenuation = distance_attenuation;

		// The point light type
		// diffuse = max(dot(N,normalize(L)),0) * attenuation
		float dotNL = (WorldNormal[0] * Lx + WorldNormal[1] * Ly + WorldNormal[2] * Lz) * rcp_dist;
		float point_attenuation = MAX(dotNL, 0.0f) * distance_attenuation;

		uint32_t attenuation = (uint32_t)(is_attenuated ? (int32_t)point_attenuation : (int32_t)simple_attenuation);

		lit[0] += (RPART(Lights[i].color) * attenuation) >> 8;
		lit[1] += (GPART(Lights[i].color) * attenuation) >> 8;
		lit[2] += (BPART(Lights[i].color) * attenuation) >> 8;
	}

	lit[0] = MIN(lit[0], (uint32_t)256);
	lit[1] = MIN(lit[1], (uint32_t)256);
	lit[2] = MIN(lit[2], (uint32_t)256);

	return MAKEARGB(lit[3], lit[0], lit[1], lit[2]);
}

/////////////////////////////////////////////////////////////////////////////

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::StepY()
{
	GradPosY.W += GradStepY.W;

	GradPosY.WorldX += GradStepY.WorldX;
	GradPosY.WorldY += GradStepY.WorldY;
	GradPosY.WorldZ += GradStepY.WorldZ;
	GradPosY.U += GradStepY.U;
	GradPosY.V += GradStepY.V;

	Dest += Pitch;
}

template<typename ModeT>
uint8_t ScreenBlockDrawer8<ModeT>::Blend(uint8_t *destptr, uint8_t srcpal, uint8_t srcalpha)
{
	using namespace TriScreenDrawerModes;

	if (ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_Zero)
	{
		return srcpal;
	}
	else if (ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_One)
	{
		uint32_t src = GPalette.BaseColors[srcpal];

		if (ModeT::BlendOp == STYLEOP_Add)
		{
			uint32_t dest = GPalette.BaseColors[*destptr];
			uint32_t out_r = MIN<uint32_t>(RPART(dest) + RPART(src), 255);
			uint32_t out_g = MIN<uint32_t>(GPART(dest) + GPART(src), 255);
			uint32_t out_b = MIN<uint32_t>(BPART(dest) + BPART(src), 255);
			return RGB256k.All[((out_r >> 2) << 12) | ((out_g >> 2) << 6) | (out_b >> 2)];
		}
		else if (ModeT::BlendOp == STYLEOP_RevSub)
		{
			uint32_t dest = GPalette.BaseColors[*destptr];
			uint32_t out_r = MAX<uint32_t>(RPART(dest) - RPART(src), 0);
			uint32_t out_g = MAX<uint32_t>(GPART(dest) - GPART(src), 0);
			uint32_t out_b = MAX<uint32_t>(BPART(dest) - BPART(src), 0);
			return RGB256k.All[((out_r >> 2) << 12) | ((out_g >> 2) << 6) | (out_b >> 2)];
		}
		else //if (ModeT::BlendOp == STYLEOP_Sub)
		{
			uint32_t dest = GPalette.BaseColors[*destptr];
			uint32_t out_r = MAX<uint32_t>(RPART(src) - RPART(dest), 0);
			uint32_t out_g = MAX<uint32_t>(GPART(src) - GPART(dest), 0);
			uint32_t out_b = MAX<uint32_t>(BPART(src) - BPART(dest), 0);
			return RGB256k.All[((out_r >> 2) << 12) | ((out_g >> 2) << 6) | (out_b >> 2)];
		}
	}
	else
	{
		if (ModeT::SWFlags & SWSTYLEF_SrcColorOneMinusSrcColor)
		{
			uint32_t src = GPalette.BaseColors[srcpal];
			uint32_t dest = GPalette.BaseColors[*destptr];
			uint32_t sfactor_r = RPART(src); sfactor_r += sfactor_r >> 7; // 255 -> 256
			uint32_t sfactor_g = GPART(src); sfactor_g += sfactor_g >> 7; // 255 -> 256
			uint32_t sfactor_b = BPART(src); sfactor_b += sfactor_b >> 7; // 255 -> 256
			uint32_t sfactor_a = srcalpha; sfactor_a += sfactor_a >> 7; // 255 -> 256
			uint32_t dfactor_r = 256 - sfactor_r;
			uint32_t dfactor_g = 256 - sfactor_g;
			uint32_t dfactor_b = 256 - sfactor_b;
			uint32_t out_r = (RPART(dest) * dfactor_r + RPART(src) * sfactor_r + 128) >> 8;
			uint32_t out_g = (GPART(dest) * dfactor_g + GPART(src) * sfactor_g + 128) >> 8;
			uint32_t out_b = (BPART(dest) * dfactor_b + BPART(src) * sfactor_b + 128) >> 8;
			return RGB256k.All[((out_r >> 2) << 12) | ((out_g >> 2) << 6) | (out_b >> 2)];
		}
		else
		{
			if (srcalpha == 255)
			{
				return srcpal;
			}
			else if (srcalpha == 0)
			{
				return *destptr;
			}

			uint32_t src = GPalette.BaseColors[srcpal];
			uint32_t dest = GPalette.BaseColors[*destptr];
			uint32_t sfactor = srcalpha; sfactor += sfactor >> 7; // 255 -> 256
			uint32_t dfactor = 256 - sfactor;
			uint32_t src_r = RPART(src) * sfactor;
			uint32_t src_g = GPART(src) * sfactor;
			uint32_t src_b = BPART(src) * sfactor;
			uint32_t dest_r = RPART(dest);
			uint32_t dest_g = GPART(dest);
			uint32_t dest_b = BPART(dest);
			if (ModeT::BlendDest == STYLEALPHA_One)
			{
				dest_r <<= 8;
				dest_g <<= 8;
				dest_b <<= 8;
			}
			else
			{
				uint32_t dfactor = 256 - sfactor;
				dest_r *= dfactor;
				dest_g *= dfactor;
				dest_b *= dfactor;
			}
			uint32_t out_r, out_g, out_b;
			if (ModeT::BlendOp == STYLEOP_Add)
			{
				if (ModeT::BlendDest == STYLEALPHA_One)
				{
					out_r = MIN<int32_t>((dest_r + src_r + 128) >> 8, 255);
					out_g = MIN<int32_t>((dest_g + src_g + 128) >> 8, 255);
					out_b = MIN<int32_t>((dest_b + src_b + 128) >> 8, 255);
				}
				else
				{
					out_r = (dest_r + src_r + 128) >> 8;
					out_g = (dest_g + src_g + 128) >> 8;
					out_b = (dest_b + src_b + 128) >> 8;
				}
			}
			else if (ModeT::BlendOp == STYLEOP_RevSub)
			{
				out_r = MAX<int32_t>(static_cast<int32_t>(dest_r - src_r + 128) >> 8, 0);
				out_g = MAX<int32_t>(static_cast<int32_t>(dest_g - src_g + 128) >> 8, 0);
				out_b = MAX<int32_t>(static_cast<int32_t>(dest_b - src_b + 128) >> 8, 0);
			}
			else //if (ModeT::BlendOp == STYLEOP_Sub)
			{
				out_r = MAX<int32_t>(static_cast<int32_t>(src_r - dest_r + 128) >> 8, 0);
				out_g = MAX<int32_t>(static_cast<int32_t>(src_g - dest_g + 128) >> 8, 0);
				out_b = MAX<int32_t>(static_cast<int32_t>(src_b - dest_b + 128) >> 8, 0);
			}
			return RGB256k.All[((out_r >> 2) << 12) | ((out_g >> 2) << 6) | (out_b >> 2)];
		}
	}
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::StoreFull(int offset)
{
	uint8_t *d = Dest + offset;
	for (int i = 0; i < 4; i++)
		d[i] = Blend(d + i, Shader.FragColor[i], Shader.FragAlpha[i]);
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::StoreMasked(int offset, uint32_t mask)
{
	using namespace TriScreenDrawerModes;

	uint8_t *d = Dest + offset;
	if (!(ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_Zero))
	{
		uint32_t m = mask;
		for (int i = 0; i < 4; i++)
		{
			if (m & (1 << 31))
				Shader.FragColor[i] = Blend(d + i, Shader.FragColor[i], Shader.FragAlpha[i]);
			m <<= 1;
		}
	}

	for (int i = 0; i < 4; i++)
	{
		if (mask & (1 << 31))
			d[i] = Shader.FragColor[i];
		mask <<= 1;
	}
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::ProcessMaskRange(uint32_t mask)
{
	for (int yy = 0; yy < 4; yy++)
	{
		GradPosX = GradPosY;

		Shader.SetVaryings(GradPosX, GradStepX);
		Shader.Run();
		StoreMasked(0, mask);
		mask <<= 4;

		Shader.SetVaryings(GradPosX, GradStepX);
		Shader.Run();
		StoreMasked(4, mask);
		mask <<= 4;

		StepY();
	}
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::ProcessBlock(uint32_t mask0, uint32_t mask1)
{
	if (mask0 == 0xffffffff && mask1 == 0xffffffff)
	{
		for (int yy = 0; yy < 8; yy++)
		{
			GradPosX = GradPosY;

			Shader.SetVaryings(GradPosX, GradStepX);
			Shader.Run();
			StoreFull(0);

			Shader.SetVaryings(GradPosX, GradStepX);
			Shader.Run();
			StoreFull(4);

			StepY();
		}
	}
	else
	{
		ProcessMaskRange(mask0);
		ProcessMaskRange(mask1);
	}
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::SetUniforms(const TriDrawTriangleArgs *args)
{
	uint32_t maskvalue = args->uniforms->FixedLight() ? 0 : 0xffffffff;
	float *maskvaluef = (float*)&maskvalue;

	Shader.LightMask = *maskvaluef;
	Shader.Light = args->uniforms->Light() * 256.0f / 255.0f;
	Shader.Shade = 2.0f - (Shader.Light + 12.0f) / 128.0f;
	Shader.GlobVis = args->uniforms->GlobVis() * (1.0f / 32.0f);
	Shader.Light /= 256.0f;
	Shader.Colormaps = args->uniforms->BaseColormap();

	Shader.Tex.data = args->uniforms->TexturePixels();
	Shader.Tex.width = args->uniforms->TextureWidth();
	Shader.Tex.height = args->uniforms->TextureHeight();

	Shader.TranslatedTex.data = args->uniforms->TexturePixels();
	Shader.TranslatedTex.translation = args->uniforms->Translation();
	Shader.TranslatedTex.width = args->uniforms->TextureWidth();
	Shader.TranslatedTex.height = args->uniforms->TextureHeight();

	Shader.SkycapColor = args->uniforms->Color();
	Shader.FillColor = args->uniforms->Color();
	Shader.Alpha = args->uniforms->SrcAlpha();

	Shader.Lights = args->uniforms->Lights();
	Shader.NumLights = args->uniforms->NumLights();
	Shader.WorldNormal[0] = args->uniforms->Normal().X;
	Shader.WorldNormal[1] = args->uniforms->Normal().Y;
	Shader.WorldNormal[2] = args->uniforms->Normal().Z;
	Shader.DynLightColor[0] = RPART(args->uniforms->DynLightColor());
	Shader.DynLightColor[1] = GPART(args->uniforms->DynLightColor());
	Shader.DynLightColor[2] = BPART(args->uniforms->DynLightColor());
	Shader.DynLightColor[3] = APART(args->uniforms->DynLightColor());
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
{
	GradStepX = gradientX;
	GradStepY = gradientY;
	GradPosY.W = v1.w + GradStepX.W * (destX - v1.x) + GradStepY.W * (destY - v1.y);
	GradPosY.U = v1.u * v1.w + GradStepX.U * (destX - v1.x) + GradStepY.U * (destY - v1.y);
	GradPosY.V = v1.v * v1.w + GradStepX.V * (destX - v1.x) + GradStepY.V * (destY - v1.y);
	GradPosY.WorldX = v1.worldX * v1.w + GradStepX.WorldX * (destX - v1.x) + GradStepY.WorldX * (destY - v1.y);
	GradPosY.WorldY = v1.worldY * v1.w + GradStepX.WorldY * (destX - v1.x) + GradStepY.WorldY * (destY - v1.y);
	GradPosY.WorldZ = v1.worldZ * v1.w + GradStepX.WorldZ * (destX - v1.x) + GradStepY.WorldZ * (destY - v1.y);
}

template<typename ModeT>
void ScreenBlockDrawer8<ModeT>::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawer8 block;
	block.Dest = args->dest + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
