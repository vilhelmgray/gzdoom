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

class TranslatedSamplerAVX2
{
public:
	uint32_t height;
	const uint8_t *data;
	const uint32_t *translation;

	__m256i size;

	FORCEINLINE void VECTORCALL SetSource(const uint8_t *pixels, const uint32_t *translation, int width, int height);
	FORCEINLINE __m256i VECTORCALL TextureNearest(const SWVec8fAVX2 &input) const;
};

template<typename ModeT>
class SWFragmentShaderAVX2
{
public:
	// Uniforms
	SamplerAVX2 Tex;
	TranslatedSamplerAVX2 TranslatedTex;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;
	PolyLight *Lights;
	int NumLights;
	__m256 WorldNormal;
	__m256i DynLightColor;
	__m256i SkycapColor;
	__m256i FillColor;
	__m256i Alpha;

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

template<typename ModeT>
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
	FORCEINLINE __m256i VECTORCALL Blend(uint32_t *destptr, __m256i src);

	// Gradients
	ScreenTriangleStepVariables GradPosX;
	ScreenTriangleStepVariables GradPosY;
	ScreenTriangleStepVariables GradStepX;
	ScreenTriangleStepVariables GradStepY;

	// Blend stage
	uint32_t *Dest;
	int Pitch;

	SWFragmentShaderAVX2<ModeT> Shader;
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

void TranslatedSamplerAVX2::SetSource(const uint8_t *pixels, const uint32_t *trans, int w, int h)
{
	data = pixels;
	translation = trans;
	height = h;
	size = _mm256_slli_epi16(_mm256_setr_epi16(w, w, w, w, h, h, h, h, w, w, w, w, h, h, h, h), 1);
}

__m256i TranslatedSamplerAVX2::TextureNearest(const SWVec8fAVX2 &input) const
{
	__m256i tmpx = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.x, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmpy = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.y, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmp = _mm256_packs_epi32(tmpx, tmpy);

	__m256i xy = _mm256_mulhi_epi16(tmp, size);
	__m256i offsetx = _mm256_mullo_epi32(_mm256_unpacklo_epi16(xy, _mm256_setzero_si256()), _mm256_set1_epi32(height));
	__m256i offsety = _mm256_unpackhi_epi16(xy, _mm256_setzero_si256());
	__m256i offset = _mm256_add_epi32(offsetx, offsety);

	uint32_t offsets[8];
	_mm256_storeu_si256((__m256i*)offsets, offset);

	return _mm256_setr_epi32(
		translation[data[offsets[0]]], translation[data[offsets[1]]], translation[data[offsets[2]]], translation[data[offsets[3]]],
		translation[data[offsets[4]]], translation[data[offsets[5]]], translation[data[offsets[6]]], translation[data[offsets[7]]]);
}

/////////////////////////////////////////////////////////////////////////////

template<typename ModeT>
void SWFragmentShaderAVX2<ModeT>::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
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

template<typename ModeT>
void SWFragmentShaderAVX2<ModeT>::Run()
{
	using namespace TriScreenDrawerModes;

	__m256i fg;
	
	if (ModeT::SWFlags & SWSTYLEF_Fill)
	{
		fg = FillColor;
	}
	else if (ModeT::SWFlags & SWSTYLEF_FogBoundary)
	{
		fg = _mm256_loadu_si256((const __m256i*)FragColor);
	}
	else if (ModeT::SWFlags & SWSTYLEF_Translated)
	{
		fg = TranslatedTex.TextureNearest(TexCoord);
	}
	else
	{
		fg = Tex.TextureNearest(TexCoord);
	}

	if ((ModeT::Flags & STYLEF_ColorIsFixed) && !(ModeT::SWFlags & SWSTYLEF_Fill))
	{
		__m256i rgbmask = _mm256_set1_epi32(0x00ffffff);
		if (ModeT::Flags & STYLEF_RedIsAlpha)
			fg = _mm256_or_si256(_mm256_andnot_si256(rgbmask, _mm256_slli_epi32(fg, 8)), _mm256_and_si256(rgbmask, FillColor));
		else
			fg = _mm256_or_si256(_mm256_andnot_si256(rgbmask, fg), _mm256_and_si256(rgbmask, FillColor));
	}

	if (!(ModeT::Flags & STYLEF_Alpha1))
	{
		__m256i a = _mm256_srli_epi32(fg, 24);
		a = _mm256_srli_epi32(_mm256_mullo_epi16(a, Alpha), 8);
		fg = _mm256_or_si256(_mm256_and_si256(fg, _mm256_set1_epi32(0x00ffffff)), _mm256_slli_epi32(a, 24));
	}

	if (ModeT::SWFlags & SWSTYLEF_Skycap)
	{
		__m256i v = _mm256_cvtps_epi32(_mm256_mul_ps(TexCoord.y, _mm256_set1_ps(1 << 24)));

		int start_fade = 2; // How fast it should fade out

		__m256i alpha_top = _mm256_max_epi32(_mm256_min_epi32(_mm256_srai_epi32(v, 16 - start_fade), _mm256_set1_epi32(256)), _mm256_setzero_si256());
		__m256i alpha_bottom = _mm256_max_epi32(_mm256_min_epi32(_mm256_srai_epi32(_mm256_sub_epi32(_mm256_set1_epi32(2 << 24), v), 16 - start_fade), _mm256_set1_epi32(256)), _mm256_setzero_si256());
		__m256i a = _mm256_min_epi32(alpha_top, alpha_bottom);
		__m256i inv_a = _mm256_sub_epi32(_mm256_set1_epi32(256), a);

		a = _mm256_or_si256(a, _mm256_slli_epi32(a, 16));
		inv_a = _mm256_or_si256(inv_a, _mm256_slli_epi32(inv_a, 16));
		__m256i alo = _mm256_shuffle_epi32(a, _MM_SHUFFLE(1, 1, 0, 0));
		__m256i ahi = _mm256_shuffle_epi32(a, _MM_SHUFFLE(3, 3, 2, 2));
		__m256i inv_alo = _mm256_shuffle_epi32(inv_a, _MM_SHUFFLE(1, 1, 0, 0));
		__m256i inv_ahi = _mm256_shuffle_epi32(inv_a, _MM_SHUFFLE(3, 3, 2, 2));

		__m256i fglo = _mm256_unpacklo_epi8(fg, _mm256_setzero_si256());
		__m256i fghi = _mm256_unpackhi_epi8(fg, _mm256_setzero_si256());

		fglo = _mm256_srli_epi16(_mm256_add_epi16(_mm256_add_epi16(_mm256_mullo_epi16(fglo, alo), _mm256_mullo_epi16(SkycapColor, inv_alo)), _mm256_set1_epi16(127)), 8);
		fghi = _mm256_srli_epi16(_mm256_add_epi16(_mm256_add_epi16(_mm256_mullo_epi16(fghi, ahi), _mm256_mullo_epi16(SkycapColor, inv_ahi)), _mm256_set1_epi16(127)), 8);

		fg = _mm256_packus_epi16(fglo, fghi);
	}
	else
	{
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

		__m256i keepalpha = _mm256_setr_epi32(0, 0, 0, 0xffffffff, 0, 0, 0, 0xffffffff);

		fg0 = _mm256_or_si256(_mm256_andnot_si256(keepalpha, _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg0), lightpos0))), _mm256_and_si256(keepalpha, fg0));
		fg1 = _mm256_or_si256(_mm256_andnot_si256(keepalpha, _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg1), lightpos1))), _mm256_and_si256(keepalpha, fg1));
		fg2 = _mm256_or_si256(_mm256_andnot_si256(keepalpha, _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg2), lightpos2))), _mm256_and_si256(keepalpha, fg2));
		fg3 = _mm256_or_si256(_mm256_andnot_si256(keepalpha, _mm256_cvtps_epi32(_mm256_mul_ps(_mm256_cvtepi32_ps(fg3), lightpos3))), _mm256_and_si256(keepalpha, fg3));

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
	}

	_mm256_storeu_si256((__m256i*)FragColor, fg);
}

template<typename ModeT>
SWVec2usAVX2 SWFragmentShaderAVX2<ModeT>::CalcDynamicLight()
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

	// Keep alpha intact
	lit.x = _mm256_or_si256(lit.x, _mm256_setr_epi16(0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255));
	lit.y = _mm256_or_si256(lit.y, _mm256_setr_epi16(0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255, 0, 0, 0, 255));

	return lit;
}

/////////////////////////////////////////////////////////////////////////////

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::StepY()
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
__m256i ScreenBlockDrawerAVX2<ModeT>::Blend(uint32_t *destptr, __m256i src)
{
	using namespace TriScreenDrawerModes;

	if (ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_Zero)
	{
		return src;
	}
	else if (ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_One)
	{
		if (ModeT::BlendOp == STYLEOP_Add)
		{
			__m256i dest = _mm256_loadu_si256((__m256i*)destptr);
			__m256i src = _mm256_loadu_si256((const __m256i*)Shader.FragColor);
			__m256i out = _mm256_adds_epu8(dest, src);
			return out;
		}
		else if (ModeT::BlendOp == STYLEOP_RevSub)
		{
			__m256i dest = _mm256_loadu_si256((__m256i*)destptr);
			__m256i src = _mm256_loadu_si256((const __m256i*)Shader.FragColor);
			__m256i out = _mm256_subs_epu8(dest, src);
			return out;
		}
		else //if (ModeT::BlendOp == STYLEOP_Sub)
		{
			__m256i dest = _mm256_loadu_si256((__m256i*)destptr);
			__m256i src = _mm256_loadu_si256((const __m256i*)Shader.FragColor);
			__m256i out = _mm256_subs_epu8(src, dest);
			return out;
		}
	}
	else
	{
		__m256i dest = _mm256_loadu_si256((__m256i*)destptr);
		__m256i destlo = _mm256_unpacklo_epi8(dest, _mm256_setzero_si256());
		__m256i desthi = _mm256_unpackhi_epi8(dest, _mm256_setzero_si256());

		__m256i srclo = _mm256_unpacklo_epi8(src, _mm256_setzero_si256());
		__m256i srchi = _mm256_unpackhi_epi8(src, _mm256_setzero_si256());

		__m256i sfactorlo, sfactorhi;
		if (ModeT::SWFlags & SWSTYLEF_SrcColorOneMinusSrcColor)
		{
			sfactorlo = srclo;
			sfactorhi = srchi;
		}
		else // if (ModeT::BlendSrc == STYLEALPHA_Src && ModeT::BlendDest == STYLEALPHA_InvSrc)
		{
			sfactorlo = _mm256_shufflehi_epi16(_mm256_shufflelo_epi16(srclo, _MM_SHUFFLE(3, 3, 3, 3)), _MM_SHUFFLE(3, 3, 3, 3));
			sfactorhi = _mm256_shufflehi_epi16(_mm256_shufflelo_epi16(srchi, _MM_SHUFFLE(3, 3, 3, 3)), _MM_SHUFFLE(3, 3, 3, 3));
		}
		sfactorlo = _mm256_add_epi16(sfactorlo, _mm256_srli_epi16(sfactorlo, 7)); // 255 -> 256
		sfactorhi = _mm256_add_epi16(sfactorhi, _mm256_srli_epi16(sfactorhi, 7)); // 255 -> 256
		srclo = _mm256_mullo_epi16(srclo, sfactorlo);
		srchi = _mm256_mullo_epi16(srchi, sfactorhi);

		if (ModeT::BlendDest == STYLEALPHA_One)
		{
			srclo = _mm256_srli_epi16(srclo, 1);
			srchi = _mm256_srli_epi16(srchi, 1);
			destlo = _mm256_slli_epi16(destlo, 7);
			desthi = _mm256_slli_epi16(desthi, 7);
		}
		else
		{
			__m256i dfactorlo = _mm256_sub_epi16(_mm256_set1_epi16(256), sfactorlo);
			__m256i dfactorhi = _mm256_sub_epi16(_mm256_set1_epi16(256), sfactorhi);
			srclo = _mm256_srli_epi16(srclo, 1);
			srchi = _mm256_srli_epi16(srchi, 1);
			destlo = _mm256_srli_epi16(_mm256_mullo_epi16(destlo, dfactorlo), 1);
			desthi = _mm256_srli_epi16(_mm256_mullo_epi16(desthi, dfactorhi), 1);
		}

		__m256i outlo, outhi;
		if (ModeT::BlendOp == STYLEOP_Add)
		{
			outlo = _mm256_adds_epi16(destlo, srclo);
			outhi = _mm256_adds_epi16(desthi, srchi);
		}
		else if (ModeT::BlendOp == STYLEOP_RevSub)
		{
			outlo = _mm256_subs_epi16(destlo, srclo);
			outhi = _mm256_subs_epi16(desthi, srchi);
		}
		else //if (ModeT::BlendOp == STYLEOP_Sub)
		{
			outlo = _mm256_subs_epi16(srclo, destlo);
			outhi = _mm256_subs_epi16(srchi, desthi);
		}

		outlo = _mm256_srai_epi16(_mm256_adds_epi16(outlo, _mm256_set1_epi16(64)), 7);
		outhi = _mm256_srai_epi16(_mm256_adds_epi16(outhi, _mm256_set1_epi16(64)), 7);
		return _mm256_packus_epi16(outlo, outhi);
	}
}

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::StoreFull()
{
	_mm256_storeu_si256((__m256i*)Dest, Blend(Dest, _mm256_loadu_si256((const __m256i*)Shader.FragColor)));
}

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::StoreMasked(uint32_t mask)
{
	using namespace TriScreenDrawerModes;

	uint32_t *d = Dest;

	if (!(ModeT::BlendSrc == STYLEALPHA_One && ModeT::BlendDest == STYLEALPHA_Zero))
	{
		uint32_t tmp[8];
		uint32_t m = mask;
		for (int i = 0; i < 8; i++)
		{
			if (m & (1 << 31))
				tmp[i] = d[i];
			m <<= 1;
		}

		_mm256_storeu_si256((__m256i*)Shader.FragColor, Blend(tmp, _mm256_loadu_si256((const __m256i*)Shader.FragColor)));
	}

	for (int i = 0; i < 8; i++)
	{
		if (mask & (1 << 31))
			d[i] = Shader.FragColor[i];
		mask <<= 1;
	}
}

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::ProcessMaskRange(uint32_t mask)
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

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::ProcessBlock(uint32_t mask0, uint32_t mask1)
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

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::SetUniforms(const TriDrawTriangleArgs *args)
{
	uint32_t maskvalue = args->uniforms->FixedLight() ? 0 : 0xffffffff;
	float *maskvaluef = (float*)&maskvalue;

	Shader.LightMask = *maskvaluef;
	Shader.Light = args->uniforms->Light() * 256.0f / 255.0f;
	Shader.Shade = 2.0f - (Shader.Light + 12.0f) / 128.0f;
	Shader.GlobVis = args->uniforms->GlobVis() * (1.0f / 32.0f);
	Shader.Light /= 256.0f;

	Shader.Tex.SetSource((const uint32_t *)args->uniforms->TexturePixels(), args->uniforms->TextureWidth(), args->uniforms->TextureHeight());
	Shader.TranslatedTex.SetSource(args->uniforms->TexturePixels(), (const uint32_t *)args->uniforms->Translation(), args->uniforms->TextureWidth(), args->uniforms->TextureHeight());

	Shader.FillColor = _mm256_set1_epi32(args->uniforms->Color());
	Shader.SkycapColor = _mm256_unpacklo_epi8(Shader.FillColor, _mm256_setzero_si256());
	Shader.Alpha = _mm256_set1_epi32(args->uniforms->SrcAlpha());

	Shader.Lights = args->uniforms->Lights();
	Shader.NumLights = args->uniforms->NumLights();
	__m128 worldNormal = _mm_setr_ps(args->uniforms->Normal().X, args->uniforms->Normal().Y, args->uniforms->Normal().Z, 0.0f);
	__m128i dynLightColor = _mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(args->uniforms->DynLightColor()), _mm_setzero_si128()), _mm_setzero_si128());
	Shader.WorldNormal = _mm256_set_m128(worldNormal, worldNormal);
	Shader.DynLightColor = _mm256_set_m128i(dynLightColor, dynLightColor);
}

template<typename ModeT>
void ScreenBlockDrawerAVX2<ModeT>::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
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
void ScreenBlockDrawerAVX2<ModeT>::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawerAVX2 block;
	block.Dest = ((uint32_t *)args->dest) + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
