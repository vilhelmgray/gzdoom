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

struct SWVec4fSSE2
{
	__m128 x, y, z, w;
};

struct SWVec2usSSE2
{
	__m128i x, y;
};

class SamplerSSE2
{
public:
	uint32_t height;
	const uint32_t *data;

	__m128i size;

	FORCEINLINE void VECTORCALL SetSource(const uint32_t *pixels, int width, int height);
	FORCEINLINE __m128i VECTORCALL TextureNearest(const SWVec4fSSE2 &input) const;
};

template<typename SamplerT>
class SWFragmentShaderSSE2
{
public:
	// Uniforms
	SamplerSSE2 Tex;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;
	PolyLight *Lights;
	int NumLights;
	__m128 WorldNormal;
	__m128i DynLightColor;
	__m128i SkycapColor;

	// In variables
	__m128 GradW;
	SWVec4fSSE2 WorldPos;
	SWVec4fSSE2 TexCoord;

	// Out variables
	uint32_t FragColor[4];

	FORCEINLINE void VECTORCALL SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void VECTORCALL Run();

private:
	FORCEINLINE SWVec2usSSE2 VECTORCALL CalcDynamicLight();
};

template<typename BlendT, typename SamplerT>
class ScreenBlockDrawerSSE2
{
public:
	static void Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args);

private:
	FORCEINLINE void VECTORCALL SetUniforms(const TriDrawTriangleArgs *args);
	FORCEINLINE void VECTORCALL SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY);
	FORCEINLINE void VECTORCALL ProcessBlock(uint32_t mask0, uint32_t mask1);
	FORCEINLINE void VECTORCALL ProcessMaskRange(uint32_t mask);
	FORCEINLINE void VECTORCALL StepY();
	FORCEINLINE void VECTORCALL StoreFull(int offset);
	FORCEINLINE void VECTORCALL StoreMasked(int offset, uint32_t mask);

	// Gradients
	ScreenTriangleStepVariables GradPosX;
	ScreenTriangleStepVariables GradPosY;
	ScreenTriangleStepVariables GradStepX;
	ScreenTriangleStepVariables GradStepY;

	// Blend stage
	uint32_t *Dest;
	int Pitch;

	SWFragmentShaderSSE2<SamplerT> Shader;
};

/////////////////////////////////////////////////////////////////////////////

void SamplerSSE2::SetSource(const uint32_t *pixels, int w, int h)
{
	data = pixels;
	height = h;
	size = _mm_slli_epi16(_mm_setr_epi16(w, w, w, w, h, h, h, h), 1);
}

__m128i SamplerSSE2::TextureNearest(const SWVec4fSSE2 &input) const
{
	__m128i tmpx = _mm_srli_epi32(_mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(input.x, _mm_set1_ps(1 << 24))), 8), 17);
	__m128i tmpy = _mm_srli_epi32(_mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(input.y, _mm_set1_ps(1 << 24))), 8), 17);
	__m128i tmp = _mm_packs_epi32(tmpx, tmpy);

	__m128i xy = _mm_mulhi_epi16(tmp, size);
#if 0 // SSE 4.1
	__m128i offsetx = _mm_mullo_epi32(_mm_unpacklo_epi16(xy, _mm_setzero_si128()), _mm_set1_epi32(height));
#else // SSE 2
	__m128i offsetx = _mm_unpacklo_epi16(xy, _mm_setzero_si128());
	__m128i mheight = _mm_set1_epi32(height);
	offsetx = _mm_or_si128(
		_mm_shuffle_epi32(_mm_mul_epu32(_mm_unpacklo_epi32(offsetx, _mm_setzero_si128()), mheight), _MM_SHUFFLE(3, 1, 2, 0)),
		_mm_shuffle_epi32(_mm_mul_epu32(_mm_unpackhi_epi32(offsetx, _mm_setzero_si128()), mheight), _MM_SHUFFLE(2, 0, 3, 1)));
#endif
	__m128i offsety = _mm_unpackhi_epi16(xy, _mm_setzero_si128());
	__m128i offset = _mm_add_epi32(offsetx, offsety);

	uint32_t offsets[4];
	_mm_storeu_si128((__m128i*)offsets, offset);

	return _mm_setr_epi32(data[offsets[0]], data[offsets[1]], data[offsets[2]], data[offsets[3]]);
}

/////////////////////////////////////////////////////////////////////////////

template<typename SamplerT>
void SWFragmentShaderSSE2<SamplerT>::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
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

	GradW = wuv0;
	TexCoord.x = _mm_mul_ps(wuv1, inv_w);
	TexCoord.y = _mm_mul_ps(wuv2, inv_w);

	__m128 stepworld = _mm_loadu_ps(&step.WorldX);
	__m128 world0 = _mm_loadu_ps(&pos.WorldX);
	__m128 world1 = _mm_add_ps(world0, stepworld);
	__m128 world2 = _mm_add_ps(world1, stepworld);
	__m128 world3 = _mm_add_ps(world2, stepworld);
	_mm_storeu_ps(&pos.WorldX, _mm_add_ps(world3, stepworld));

	_MM_TRANSPOSE4_PS(world0, world1, world2, world3);

	WorldPos.x = _mm_mul_ps(world0, inv_w);
	WorldPos.y = _mm_mul_ps(world1, inv_w);
	WorldPos.z = _mm_mul_ps(world2, inv_w);
}

template<typename SamplerT>
void SWFragmentShaderSSE2<SamplerT>::Run()
{
	using namespace TriScreenDrawerModes;

	__m128i fg = Tex.TextureNearest(TexCoord);

	if (SamplerT::Mode == (int)Samplers::Skycap)
	{
		__m128i v = _mm_cvtps_epi32(_mm_mul_ps(TexCoord.y, _mm_set1_ps(1 << 24)));

		int start_fade = 2; // How fast it should fade out

		__m128i alpha_top = _mm_max_epi16(_mm_min_epi16(_mm_srai_epi32(v, 16 - start_fade), _mm_set1_epi32(256)), _mm_setzero_si128());
		__m128i alpha_bottom = _mm_max_epi16(_mm_min_epi16(_mm_srai_epi32(_mm_sub_epi32(_mm_set1_epi32(2 << 24), v), 16 - start_fade), _mm_set1_epi32(256)), _mm_setzero_si128());
		__m128i a = _mm_and_si128(_mm_min_epi16(alpha_top, alpha_bottom), _mm_set1_epi32(0xffff));
		__m128i inv_a = _mm_sub_epi32(_mm_set1_epi32(256), a);

		a = _mm_or_si128(a, _mm_slli_epi32(a, 16));
		inv_a = _mm_or_si128(inv_a, _mm_slli_epi32(inv_a, 16));
		__m128i alo = _mm_shuffle_epi32(a, _MM_SHUFFLE(1, 1, 0, 0));
		__m128i ahi = _mm_shuffle_epi32(a, _MM_SHUFFLE(3, 3, 2, 2));
		__m128i inv_alo = _mm_shuffle_epi32(inv_a, _MM_SHUFFLE(1, 1, 0, 0));
		__m128i inv_ahi = _mm_shuffle_epi32(inv_a, _MM_SHUFFLE(3, 3, 2, 2));

		__m128i fglo = _mm_unpacklo_epi8(fg, _mm_setzero_si128());
		__m128i fghi = _mm_unpackhi_epi8(fg, _mm_setzero_si128());

		fglo = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(fglo, alo), _mm_mullo_epi16(SkycapColor, inv_alo)), _mm_set1_epi16(127)), 8);
		fghi = _mm_srli_epi16(_mm_add_epi16(_mm_add_epi16(_mm_mullo_epi16(fghi, ahi), _mm_mullo_epi16(SkycapColor, inv_ahi)), _mm_set1_epi16(127)), 8);

		fg = _mm_packus_epi16(fglo, fghi);
	}
	else
	{
		__m128 mlightmask = _mm_set1_ps(LightMask);
		__m128 lightposf = _mm_sub_ps(_mm_set1_ps(Shade), _mm_min_ps(_mm_set1_ps(24.0f / 32.0f), _mm_mul_ps(_mm_set1_ps(GlobVis), GradW)));
		lightposf = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_max_ps(_mm_min_ps(lightposf, _mm_set1_ps(31.0f / 32.0f)), _mm_setzero_ps()));
		lightposf = _mm_or_ps(_mm_and_ps(mlightmask, lightposf), _mm_andnot_ps(mlightmask, _mm_set1_ps(Light)));

		__m128 lightpos0 = _mm_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(0, 0, 0, 0));
		__m128 lightpos1 = _mm_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(1, 1, 1, 1));
		__m128 lightpos2 = _mm_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(2, 2, 2, 2));
		__m128 lightpos3 = _mm_shuffle_ps(lightposf, lightposf, _MM_SHUFFLE(3, 3, 3, 3));

		__m128i fglo = _mm_unpacklo_epi8(fg, _mm_setzero_si128());
		__m128i fghi = _mm_unpackhi_epi8(fg, _mm_setzero_si128());
		__m128i fg0 = _mm_unpacklo_epi16(fglo, _mm_setzero_si128());
		__m128i fg1 = _mm_unpackhi_epi16(fglo, _mm_setzero_si128());
		__m128i fg2 = _mm_unpacklo_epi16(fghi, _mm_setzero_si128());
		__m128i fg3 = _mm_unpackhi_epi16(fghi, _mm_setzero_si128());

		fg0 = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fg0), lightpos0));
		fg1 = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fg1), lightpos1));
		fg2 = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fg2), lightpos2));
		fg3 = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fg3), lightpos3));

		if (NumLights == 0)
		{
			fglo = _mm_packs_epi32(fg0, fg1);
			fghi = _mm_packs_epi32(fg2, fg3);
			fg = _mm_packus_epi16(fglo, fghi);
		}
		else
		{
			SWVec2usSSE2 dynlight = CalcDynamicLight();
			fglo = _mm_add_epi16(_mm_packs_epi32(fg0, fg1), _mm_srli_epi16(_mm_mullo_epi16(fglo, dynlight.x), 8));
			fghi = _mm_add_epi16(_mm_packs_epi32(fg2, fg3), _mm_srli_epi16(_mm_mullo_epi16(fghi, dynlight.y), 8));
			fg = _mm_packus_epi16(fglo, fghi);
		}
	}

	_mm_storeu_si128((__m128i*)FragColor, fg);
}

template<typename SamplerT>
SWVec2usSSE2 SWFragmentShaderSSE2<SamplerT>::CalcDynamicLight()
{
	SWVec2usSSE2 lit;

	lit.x = DynLightColor;
	lit.y = DynLightColor;

	__m128 m256f = _mm_set1_ps(256.0f);
	__m128i m256i = _mm_set1_epi16(256);
	__m128 mSignBit = _mm_set1_ps(-0.0f);

	__m128 worldnormalx = _mm_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(0, 0, 0, 0));
	__m128 worldnormaly = _mm_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(1, 1, 1, 1));
	__m128 worldnormalz = _mm_shuffle_ps(WorldNormal, WorldNormal, _MM_SHUFFLE(2, 2, 2, 2));

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

		// L = light-pos
		// dist = sqrt(dot(L, L))
		// distance_attenuation = 1 - MIN(dist * (1/radius), 1)
		__m128 Lx = _mm_sub_ps(lightposx, WorldPos.x);
		__m128 Ly = _mm_sub_ps(lightposy, WorldPos.y);
		__m128 Lz = _mm_sub_ps(lightposz, WorldPos.z);
		__m128 dist2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(Lx, Lx), _mm_mul_ps(Ly, Ly)), _mm_mul_ps(Lz, Lz));
		__m128 rcp_dist = _mm_rsqrt_ps(dist2);
		__m128 dist = _mm_mul_ps(dist2, rcp_dist);
		__m128 distance_attenuation = _mm_sub_ps(m256f, _mm_min_ps(_mm_mul_ps(dist, light_radius), m256f));

		// The simple light type
		__m128 simple_attenuation = distance_attenuation;

		// The point light type
		// diffuse = max(dot(N,normalize(L)),0) * attenuation
		__m128 dotNL = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(worldnormalx, Lx), _mm_mul_ps(worldnormaly, Ly)), _mm_mul_ps(worldnormalz, Lz)), rcp_dist);
		__m128 point_attenuation = _mm_mul_ps(_mm_max_ps(dotNL, _mm_setzero_ps()), distance_attenuation);

		__m128i attenuation = _mm_cvtps_epi32(_mm_or_ps(_mm_and_ps(is_attenuated, simple_attenuation), _mm_andnot_ps(is_attenuated, point_attenuation)));

		__m128i light_color = _mm_cvtsi32_si128(Lights[i].color);
		light_color = _mm_unpacklo_epi8(light_color, _mm_setzero_si128());
		light_color = _mm_shuffle_epi32(light_color, _MM_SHUFFLE(1, 0, 1, 0));

		__m128i attenuationlo = _mm_packs_epi32(_mm_shuffle_epi32(attenuation, _MM_SHUFFLE(0, 0, 0, 0)), _mm_shuffle_epi32(attenuation, _MM_SHUFFLE(1, 1, 1, 1)));
		__m128i attenuationhi = _mm_packs_epi32(_mm_shuffle_epi32(attenuation, _MM_SHUFFLE(2, 2, 2, 2)), _mm_shuffle_epi32(attenuation, _MM_SHUFFLE(3, 3, 3, 3)));

		lit.x = _mm_add_epi16(_mm_srli_epi16(_mm_mullo_epi16(light_color, attenuationlo), 8), lit.x);
		lit.y = _mm_add_epi16(_mm_srli_epi16(_mm_mullo_epi16(light_color, attenuationhi), 8), lit.y);
	}

	lit.x = _mm_min_epi16(lit.x, m256i);
	lit.y = _mm_min_epi16(lit.y, m256i);

	return lit;
}

/////////////////////////////////////////////////////////////////////////////

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::StepY()
{
	GradPosY.W += GradStepY.W;

	GradPosY.WorldX += GradStepY.WorldX;
	GradPosY.WorldY += GradStepY.WorldY;
	GradPosY.WorldZ += GradStepY.WorldZ;
	GradPosY.U += GradStepY.U;
	GradPosY.V += GradStepY.V;

	Dest += Pitch;
}

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::StoreFull(int offset)
{
	using namespace TriScreenDrawerModes;

	if (BlendT::Mode == (int)BlendModes::Opaque)
	{
		_mm_storeu_si128((__m128i*)(Dest + offset), _mm_loadu_si128((const __m128i*)Shader.FragColor));
	}
	else if (BlendT::Mode == (int)BlendModes::Masked)
	{
		__m128i fragcolor = _mm_loadu_si128((const __m128i*)Shader.FragColor);
		__m128i mask = _mm_cmpeq_epi32(_mm_and_si128(fragcolor, _mm_set1_epi32(0xff000000)), _mm_setzero_si128());
		mask = _mm_xor_si128(mask, _mm_set1_epi32(0xffffffff));
#if 0 // AVX 2
		_mm_maskstore_epi32((int*)Dest, mask, fragcolor);
#else // SSE 2
		_mm_maskmoveu_si128(fragcolor, mask, (char*)(Dest + offset));
#endif
	}
}

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::StoreMasked(int offset, uint32_t mask)
{
	using namespace TriScreenDrawerModes;

	if (BlendT::Mode == (int)BlendModes::Opaque)
	{
		uint32_t *d = Dest + offset;
		for (int i = 0; i < 4; i++)
		{
			if (mask & (1 << 31))
				d[i] = Shader.FragColor[i];
			mask <<= 1;
		}
	}
	else if (BlendT::Mode == (int)BlendModes::Masked)
	{
		uint32_t *d = Dest + offset;
		for (int i = 0; i < 4; i++)
		{
			if ((mask & (1 << 31)) && APART(Shader.FragColor[i]))
				d[i] = Shader.FragColor[i];
			mask <<= 1;
		}
	}
}

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::ProcessMaskRange(uint32_t mask)
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

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::ProcessBlock(uint32_t mask0, uint32_t mask1)
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

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::SetUniforms(const TriDrawTriangleArgs *args)
{
	uint32_t maskvalue = args->uniforms->FixedLight() ? 0 : 0xffffffff;
	float *maskvaluef = (float*)&maskvalue;

	Shader.LightMask = *maskvaluef;
	Shader.Light = args->uniforms->Light() * 256.0f / 255.0f;
	Shader.Shade = 2.0f - (Shader.Light + 12.0f) / 128.0f;
	Shader.GlobVis = args->uniforms->GlobVis() * (1.0f / 32.0f);
	Shader.Light /= 256.0f;

	Shader.Tex.SetSource((const uint32_t *)args->uniforms->TexturePixels(), args->uniforms->TextureWidth(), args->uniforms->TextureHeight());

	Shader.SkycapColor = _mm_unpacklo_epi8(_mm_set1_epi32(args->uniforms->Color()), _mm_setzero_si128());

	Shader.Lights = args->uniforms->Lights();
	Shader.NumLights = args->uniforms->NumLights();
	Shader.WorldNormal = _mm_setr_ps(args->uniforms->Normal().X, args->uniforms->Normal().Y, args->uniforms->Normal().Z, 0.0f);
	Shader.DynLightColor = _mm_unpacklo_epi16(_mm_unpacklo_epi8(_mm_cvtsi32_si128(args->uniforms->DynLightColor()), _mm_setzero_si128()), _mm_setzero_si128());
}

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
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

template<typename BlendT, typename SamplerT>
void ScreenBlockDrawerSSE2<BlendT, SamplerT>::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawerSSE2 block;
	block.Dest = ((uint32_t *)args->dest) + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
