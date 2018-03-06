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
#include <immintrin.h>

struct SWVec8fAVX2
{
	__m256 x, y, z, w;
};

struct SWVec2usAVX2
{
	__m256i x, y;
};

class SamplerAVX2
{
public:
	uint32_t height;
	const uint32_t *data;

	__m256i size;

	FORCEINLINE void VECTORCALL SetSource(const uint32_t *pixels, int width, int height);
	FORCEINLINE __m256i VECTORCALL TextureNearest(const SWVec8fAVX2 &input) const;
};

class SWFragmentShaderAVX2
{
public:
	// Uniforms
	SamplerAVX2 Tex;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;
	PolyLight *Lights;
	int NumLights;
	__m256 WorldNormal;
	__m256i DynLightColor;

	// In variables
	__m256 GradW;
	SWVec8fAVX2 WorldPos;
	SWVec8fAVX2 TexCoord;

	// Out variables
	uint32_t FragColor[8];

	FORCEINLINE void VECTORCALL SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void VECTORCALL Run();

private:
	FORCEINLINE SWVec2usAVX2 VECTORCALL CalcDynamicLight();
};

class ScreenBlockDrawerAVX2
{
public:
	static void Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args);

private:
	FORCEINLINE void VECTORCALL SetUniforms(const TriDrawTriangleArgs *args);
	FORCEINLINE void VECTORCALL SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY);
	FORCEINLINE void VECTORCALL ProcessBlock(uint32_t mask0, uint32_t mask1);
	FORCEINLINE void VECTORCALL ProcessMaskRange(uint32_t mask);
	FORCEINLINE void VECTORCALL StepY();
	FORCEINLINE void VECTORCALL StoreFull();
	FORCEINLINE void VECTORCALL StoreMasked(uint32_t mask);

	// Gradients
	ScreenTriangleStepVariables GradPosX;
	ScreenTriangleStepVariables GradPosY;
	ScreenTriangleStepVariables GradStepX;
	ScreenTriangleStepVariables GradStepY;

	// Blend stage
	uint32_t *Dest;
	int Pitch;

	SWFragmentShaderAVX2 Shader;
};

/////////////////////////////////////////////////////////////////////////////

void SamplerAVX2::SetSource(const uint32_t *pixels, int w, int h)
{
	data = pixels;
	height = h;
	size = _mm256_slli_epi16(_mm256_setr_epi16(w, w, w, w, h, h, h, h, w, w, w, w, h, h, h, h), 1);
}

__m256i SamplerAVX2::TextureNearest(const SWVec8fAVX2 &input) const
{
	__m256i tmpx = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.x, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmpy = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.y, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmp = _mm256_packs_epi32(tmpx, tmpy);

	__m256i xy = _mm256_mulhi_epi16(tmp, size);
	__m256i offsetx = _mm256_mullo_epi32(_mm256_unpacklo_epi16(xy, _mm256_setzero_si256()), _mm256_set1_epi32(height));
	__m256i offsety = _mm256_unpackhi_epi16(xy, _mm256_setzero_si256());
	__m256i offset = _mm256_add_epi32(offsetx, offsety);

	return _mm256_i32gather_epi32((const int*)data, offset, 4);
}

/////////////////////////////////////////////////////////////////////////////

void SWFragmentShaderAVX2::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
{
	{
		__m128 stepwuv = _mm_loadu_ps(&step.W);
		__m128 wuv0 = _mm_loadu_ps(&pos.W);
		__m128 wuv1 = _mm_add_ps(wuv0, stepwuv);
		__m128 wuv2 = _mm_add_ps(wuv1, stepwuv);
		__m128 wuv3 = _mm_add_ps(wuv2, stepwuv);
		_mm_storeu_ps(&pos.W, _mm_add_ps(wuv3, stepwuv));

		_MM_TRANSPOSE4_PS(wuv0, wuv1, wuv2, wuv3);

		//__m128 inv_w = _mm_rcp_ps(wuv0); // Much faster, but also lower quality
		__m128 inv_w = _mm_div_ps(_mm_set1_ps(1.0f), wuv0);

		GradW = _mm256_insertf128_ps(GradW, wuv0, 0);
		TexCoord.x = _mm256_insertf128_ps(TexCoord.x, _mm_mul_ps(wuv1, inv_w), 0);
		TexCoord.y = _mm256_insertf128_ps(TexCoord.y, _mm_mul_ps(wuv2, inv_w), 0);

		__m128 stepworld = _mm_loadu_ps(&step.WorldX);
		__m128 world0 = _mm_loadu_ps(&pos.WorldX);
		__m128 world1 = _mm_add_ps(world0, stepworld);
		__m128 world2 = _mm_add_ps(world1, stepworld);
		__m128 world3 = _mm_add_ps(world2, stepworld);
		_mm_storeu_ps(&pos.WorldX, _mm_add_ps(world3, stepworld));

		_MM_TRANSPOSE4_PS(world0, world1, world2, world3);

		WorldPos.x = _mm256_insertf128_ps(WorldPos.x, _mm_mul_ps(world0, inv_w), 0);
		WorldPos.y = _mm256_insertf128_ps(WorldPos.y, _mm_mul_ps(world1, inv_w), 0);
		WorldPos.z = _mm256_insertf128_ps(WorldPos.z, _mm_mul_ps(world2, inv_w), 0);
	}

	{
		__m128 stepwuv = _mm_loadu_ps(&step.W);
		__m128 wuv0 = _mm_loadu_ps(&pos.W);
		__m128 wuv1 = _mm_add_ps(wuv0, stepwuv);
		__m128 wuv2 = _mm_add_ps(wuv1, stepwuv);
		__m128 wuv3 = _mm_add_ps(wuv2, stepwuv);
		_mm_storeu_ps(&pos.W, _mm_add_ps(wuv3, stepwuv));

		_MM_TRANSPOSE4_PS(wuv0, wuv1, wuv2, wuv3);

		//__m128 inv_w = _mm_rcp_ps(wuv0); // Much faster, but also lower quality
		__m128 inv_w = _mm_div_ps(_mm_set1_ps(1.0f), wuv0);

		GradW = _mm256_insertf128_ps(GradW, wuv0, 1);
		TexCoord.x = _mm256_insertf128_ps(TexCoord.x, _mm_mul_ps(wuv1, inv_w), 1);
		TexCoord.y = _mm256_insertf128_ps(TexCoord.y, _mm_mul_ps(wuv2, inv_w), 1);

		__m128 stepworld = _mm_loadu_ps(&step.WorldX);
		__m128 world0 = _mm_loadu_ps(&pos.WorldX);
		__m128 world1 = _mm_add_ps(world0, stepworld);
		__m128 world2 = _mm_add_ps(world1, stepworld);
		__m128 world3 = _mm_add_ps(world2, stepworld);
		_mm_storeu_ps(&pos.WorldX, _mm_add_ps(world3, stepworld));

		_MM_TRANSPOSE4_PS(world0, world1, world2, world3);

		WorldPos.x = _mm256_insertf128_ps(WorldPos.x, _mm_mul_ps(world0, inv_w), 1);
		WorldPos.y = _mm256_insertf128_ps(WorldPos.y, _mm_mul_ps(world1, inv_w), 1);
		WorldPos.z = _mm256_insertf128_ps(WorldPos.z, _mm_mul_ps(world2, inv_w), 1);
	}
}

void SWFragmentShaderAVX2::Run()
{
	__m256i fg = Tex.TextureNearest(TexCoord);

	__m256 mlightmask = _mm256_set1_ps(LightMask);
	__m256 lightposf = _mm256_sub_ps(_mm256_set1_ps(Shade), _mm256_min_ps(_mm256_set1_ps(24.0f / 32.0f), _mm256_mul_ps(_mm256_set1_ps(GlobVis), GradW)));
	lightposf = _mm256_sub_ps(_mm256_set1_ps(1.0f), _mm256_max_ps(_mm256_min_ps(lightposf, _mm256_set1_ps(31.0f / 32.0f)), _mm256_setzero_ps()));
	lightposf = _mm256_or_ps(_mm256_and_ps(mlightmask, lightposf), _mm256_andnot_ps(mlightmask, _mm256_set1_ps(Light)));

	__m256 lightpos0 = _mm256_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(0, 0, 0, 0));
	__m256 lightpos1 = _mm256_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(1, 1, 1, 1));
	__m256 lightpos2 = _mm256_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(2, 2, 2, 2));
	__m256 lightpos3 = _mm256_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(3, 3, 3, 3));

	__m256i fglo = _mm256_unpacklo_epi8(fg, _mm256_setzero_si256());
	__m256i fghi = _mm256_unpackhi_epi8(fg, _mm256_setzero_si256());
	__m256i fg0 = _mm256_unpacklo_epi16(fglo, _mm256_setzero_si256());
	__m256i fg1 = _mm256_unpackhi_epi16(fglo, _mm256_setzero_si256());
	__m256i fg2 = _mm256_unpacklo_epi16(fghi, _mm256_setzero_si256());
	__m256i fg3 = _mm256_unpackhi_epi16(fghi, _mm256_setzero_si256());

	fg0 = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg0), lightpos0));
	fg1 = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg1), lightpos1));
	fg2 = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg2), lightpos2));
	fg3 = _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg3), lightpos3));

	if (NumLights == 0)
	{
		fglo = _mm256_packs_epi32(fg0, fg1);
		fghi = _mm256_packs_epi32(fg2, fg3);
		fg = _mm256_packus_epi16(fglo, fghi);
	}
	else
	{
		SWVec2usAVX2 dynlight = CalcDynamicLight();
		fglo = _mm256_add_epi16(_mm256_packs_epi32(fg0, fg1), _mm256_srli_epi16(_mm256_mullo_epi16(fglo, dynlight.x), 8));
		fghi = _mm256_add_epi16(_mm256_packs_epi32(fg2, fg3), _mm256_srli_epi16(_mm256_mullo_epi16(fghi, dynlight.y), 8));
		fg = _mm256_packus_epi16(fglo, fghi);
	}

	_mm256_storeu_si256((__m256i*)FragColor, fg);
}

SWVec2usAVX2 SWFragmentShaderAVX2::CalcDynamicLight()
{
	SWVec2usAVX2 lit;

	lit.x = DynLightColor;
	lit.y = DynLightColor;

	__m256 m256f = _mm256_set1_ps(256.0f);
	__m256i m256i = _mm256_set1_epi16(256);
	__m128 mSignBit = _mm_set1_ps(-0.0f);

	__m256 worldnormalx = _mm256_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(0, 0, 0, 0));
	__m256 worldnormaly = _mm256_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(1, 1, 1, 1));
	__m256 worldnormalz = _mm256_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(2, 2, 2, 2));

	for (int i = 0; i < NumLights; i++)
	{
		__m128 lightpos = _mm_loadu_ps(&Lights[i].x);
		__m128 lightposx = _mm_shuffle_ps(lightpos, lightpos, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 lightposy = _mm_shuffle_ps(lightpos, lightpos, _MM_SHUFFLE(1, 1, 1, 1));
		__m128 lightposz = _mm_shuffle_ps(lightpos, lightpos, _MM_SHUFFLE(2, 2, 2, 2));
		__m128 light_radius = _mm_shuffle_ps(lightpos, lightpos, _MM_SHUFFLE(3, 3, 3, 3)); // Lights[i].radius

		__m128 is_attenuated = _mm_cmpge_ss(light_radius, _mm_setzero_ps());
		is_attenuated = _mm_shuffle_ps(is_attenuated, is_attenuated, _MM_SHUFFLE(0, 0, 0, 0));
		light_radius = _mm_andnot_ps(mSignBit, light_radius);

		__m256 is_attenuated256 = _mm256_set_m128(is_attenuated, is_attenuated);

		// L = light-pos
		// dist = sqrt(dot(L, L))
		// distance_attenuation = 1 - MIN(dist * (1/radius), 1)
		__m256 Lx = _mm256_sub_ps(_mm256_set_m128(lightposx, lightposx), WorldPos.x);
		__m256 Ly = _mm256_sub_ps(_mm256_set_m128(lightposy, lightposy), WorldPos.y);
		__m256 Lz = _mm256_sub_ps(_mm256_set_m128(lightposz, lightposz), WorldPos.z);
		__m256 dist2 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(Lx, Lx), _mm256_mul_ps(Ly, Ly)), _mm256_mul_ps(Lz, Lz));
		__m256 rcp_dist = _mm256_rsqrt_ps(dist2);
		__m256 dist = _mm256_mul_ps(dist2, rcp_dist);
		__m256 distance_attenuation = _mm256_sub_ps(m256f, _mm256_min_ps(_mm256_mul_ps(dist, _mm256_set_m128(light_radius, light_radius)), m256f));

		// The simple light type
		__m256 simple_attenuation = distance_attenuation;

		// The point light type
		// diffuse = max(dot(N,normalize(L)),0) * attenuation
		__m256 dotNL = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(worldnormalx, Lx), _mm256_mul_ps(worldnormaly, Ly)), _mm256_mul_ps(worldnormalz, Lz)), rcp_dist);
		__m256 point_attenuation = _mm256_mul_ps(_mm256_max_ps(dotNL, _mm256_setzero_ps()), distance_attenuation);

		__m256i attenuation = _mm256_cvtps_epi32(_mm256_or_ps(_mm256_and_ps(is_attenuated256, simple_attenuation), _mm256_andnot_ps(is_attenuated256, point_attenuation)));

		__m128i light_color = _mm_cvtsi32_si128(Lights[i].color);
		light_color = _mm_unpacklo_epi8(light_color, _mm_setzero_si128());
		light_color = _mm_shuffle_epi32(light_color, _MM_SHUFFLE(1, 0, 1, 0));
		__m256i light_color256 = _mm256_set_m128i(light_color, light_color);

		__m256i attenuationlo = _mm256_packs_epi32(_mm256_shuffle_epi32(attenuation, _MM_SHUFFLE(0, 0, 0, 0)), _mm256_shuffle_epi32(attenuation, _MM_SHUFFLE(1, 1, 1, 1)));
		__m256i attenuationhi = _mm256_packs_epi32(_mm256_shuffle_epi32(attenuation, _MM_SHUFFLE(2, 2, 2, 2)), _mm256_shuffle_epi32(attenuation, _MM_SHUFFLE(3, 3, 3, 3)));

		lit.x = _mm256_add_epi16(_mm256_srli_epi16(_mm256_mullo_epi16(light_color256, attenuationlo), 8), lit.x);
		lit.y = _mm256_add_epi16(_mm256_srli_epi16(_mm256_mullo_epi16(light_color256, attenuationhi), 8), lit.y);
	}

	lit.x = _mm256_min_epi16(lit.x, m256i);
	lit.y = _mm256_min_epi16(lit.y, m256i);

	return lit;
}

/////////////////////////////////////////////////////////////////////////////

void ScreenBlockDrawerAVX2::StepY()
{
	GradPosY.W += GradStepY.W;

	GradPosY.WorldX += GradStepY.WorldX;
	GradPosY.WorldY += GradStepY.WorldY;
	GradPosY.WorldZ += GradStepY.WorldZ;
	GradPosY.U += GradStepY.U;
	GradPosY.V += GradStepY.V;

	Dest += Pitch;
}

void ScreenBlockDrawerAVX2::StoreFull()
{
	_mm256_storeu_si256((__m256i*)Dest, _mm256_loadu_si256((const __m256i*)Shader.FragColor));
}

void ScreenBlockDrawerAVX2::StoreMasked(uint32_t mask)
{
	uint32_t *d = Dest;
	for (int i = 0; i < 8; i++)
	{
		if (mask & (1 << 31))
			d[i] = Shader.FragColor[i];
		mask <<= 1;
	}
}

void ScreenBlockDrawerAVX2::ProcessMaskRange(uint32_t mask)
{
	for (int yy = 0; yy < 4; yy++)
	{
		GradPosX = GradPosY;

		Shader.SetVaryings(GradPosX, GradStepX);
		Shader.Run();
		StoreMasked(mask);
		mask <<= 8;

		StepY();
	}
}

void ScreenBlockDrawerAVX2::ProcessBlock(uint32_t mask0, uint32_t mask1)
{
	if (mask0 == 0xffffffff && mask1 == 0xffffffff)
	{
		for (int yy = 0; yy < 8; yy++)
		{
			GradPosX = GradPosY;

			Shader.SetVaryings(GradPosX, GradStepX);
			Shader.Run();
			StoreFull();

			StepY();
		}
	}
	else
	{
		ProcessMaskRange(mask0);
		ProcessMaskRange(mask1);
	}
}

void ScreenBlockDrawerAVX2::SetUniforms(const TriDrawTriangleArgs *args)
{
	uint32_t maskvalue = args->uniforms->FixedLight() ? 0 : 0xffffffff;
	float *maskvaluef = (float*)&maskvalue;

	Shader.LightMask = *maskvaluef;
	Shader.Light = args->uniforms->Light() * 256.0f / 255.0f;
	Shader.Shade = 2.0f - (Shader.Light + 12.0f) / 128.0f;
	Shader.GlobVis = args->uniforms->GlobVis() * (1.0f / 32.0f);
	Shader.Light /= 256.0f;

	Shader.Tex.SetSource((const uint32_t *)args->uniforms->TexturePixels(), args->uniforms->TextureWidth(), args->uniforms->TextureHeight());

	Shader.Lights = args->uniforms->Lights();
	Shader.NumLights = args->uniforms->NumLights();
	__m128 worldNormal = _mm_setr_ps(args->uniforms->Normal().X, args->uniforms->Normal().Y, args->uniforms->Normal().Z, 0.0f);
	__m128i dynLightColor = _mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(args->uniforms->DynLightColor()), _mm_setzero_si128()), _mm_setzero_si128());
	Shader.WorldNormal = _mm256_set_m128(worldNormal, worldNormal);
	Shader.DynLightColor = _mm256_set_m128i(dynLightColor, dynLightColor);
}

void ScreenBlockDrawerAVX2::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
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

void ScreenBlockDrawerAVX2::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawerAVX2 block;
	block.Dest = ((uint32_t *)args->dest) + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
