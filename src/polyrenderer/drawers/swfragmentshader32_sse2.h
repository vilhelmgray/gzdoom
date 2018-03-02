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

class SamplerSSE2
{
public:
	uint32_t width;
	uint32_t height;
	const uint32_t *data;

	FORCEINLINE __m128i TextureNearest(const SWVec4fSSE2 &input) const;
};

class SWFragmentShaderSSE2
{
public:
	// Uniforms
	SamplerSSE2 Tex;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;

	// In variables
	__m128 GradW;
	SWVec4fSSE2 WorldPos;
	SWVec4fSSE2 TexCoord;

	// Out variables
	uint32_t FragColor[4];

	FORCEINLINE void SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void Run();
};

class ScreenBlockDrawerSSE2
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

	// Gradients
	ScreenTriangleStepVariables GradPosX;
	ScreenTriangleStepVariables GradPosY;
	ScreenTriangleStepVariables GradStepX;
	ScreenTriangleStepVariables GradStepY;

	// Blend stage
	uint32_t *Dest;
	int Pitch;

	SWFragmentShaderSSE2 Shader;
};

/////////////////////////////////////////////////////////////////////////////

__m128i SamplerSSE2::TextureNearest(const SWVec4fSSE2 &input) const
{
	__m128i tmpx = _mm_srli_epi32(_mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(input.x, _mm_set1_ps(1 << 24))), 8), 17);
	__m128i tmpy = _mm_srli_epi32(_mm_slli_epi32(_mm_cvtps_epi32(_mm_mul_ps(input.y, _mm_set1_ps(1 << 24))), 8), 17);
	__m128i tmp = _mm_packs_epi32(tmpx, tmpy);
	__m128i size = _mm_setr_epi16(width << 1, width << 1, width << 1, width << 1, height << 1, height << 1, height << 1, height << 1);
	uint16_t xy[8];
	_mm_storeu_si128((__m128i*)xy, _mm_mulhi_epi16(tmp, size));

	uint32_t pixel[4];
	pixel[0] = data[xy[4] + xy[0] * height];
	pixel[1] = data[xy[5] + xy[1] * height];
	pixel[2] = data[xy[6] + xy[2] * height];
	pixel[3] = data[xy[7] + xy[3] * height];

	return _mm_loadu_si128((const __m128i*)pixel);
}

/////////////////////////////////////////////////////////////////////////////

void SWFragmentShaderSSE2::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
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

#if 0
void SWFragmentShaderSSE2::Run()
{
	__m128i fg = Tex.TextureNearest(TexCoord);

	__m128i fgX, fgY, fgZ, fgW;
	{
		__m128i mpixello = _mm_unpacklo_epi8(fg, _mm_setzero_si128());
		__m128i mpixelhi = _mm_unpackhi_epi8(fg, _mm_setzero_si128());
		__m128 p0 = _mm_castsi128_ps(_mm_unpacklo_epi16(mpixello, _mm_setzero_si128()));
		__m128 p1 = _mm_castsi128_ps(_mm_unpackhi_epi16(mpixello, _mm_setzero_si128()));
		__m128 p2 = _mm_castsi128_ps(_mm_unpacklo_epi16(mpixelhi, _mm_setzero_si128()));
		__m128 p3 = _mm_castsi128_ps(_mm_unpackhi_epi16(mpixelhi, _mm_setzero_si128()));
		_MM_TRANSPOSE4_PS(p0, p1, p2, p3);
		fgX = _mm_castps_si128(p2);
		fgY = _mm_castps_si128(p1);
		fgZ = _mm_castps_si128(p0);
		fgW = _mm_castps_si128(p3);
	}

	__m128 lightposf = _mm_sub_ps(_mm_set1_ps(Shade), _mm_min_ps(_mm_set1_ps(24.0f / 32.0f), _mm_mul_ps(_mm_set1_ps(GlobVis), _mm_load_ps(GradW))));
	lightposf = _mm_sub_ps(_mm_set1_ps(1.0f), _mm_max_ps(_mm_min_ps(lightposf, _mm_set1_ps(31.0f / 32.0f)), _mm_setzero_ps()));

	__m128 mlightmask = _mm_set1_ps(LightMask);
	lightposf = _mm_or_ps(_mm_and_ps(mlightmask, lightposf), _mm_andnot_ps(mlightmask, _mm_set1_ps(Light)));

	__m128i x = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fgX), lightposf));
	__m128i y = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fgY), lightposf));
	__m128i z = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(fgZ), lightposf));
	__m128i w = fgW;

	__m128i ppacked;
	{
		__m128 p0 = _mm_castsi128_ps(z);
		__m128 p1 = _mm_castsi128_ps(y);
		__m128 p2 = _mm_castsi128_ps(x);
		__m128 p3 = _mm_castsi128_ps(w);
		_MM_TRANSPOSE4_PS(p0, p1, p2, p3);
		ppacked = _mm_packus_epi16(_mm_packs_epi32(_mm_castps_si128(p0), _mm_castps_si128(p1)), _mm_packs_epi32(_mm_castps_si128(p2), _mm_castps_si128(p3)));
	}
	_mm_storeu_si128((__m128i*)FragColor, ppacked);
}
#else
void SWFragmentShaderSSE2::Run()
{
	__m128i fg = Tex.TextureNearest(TexCoord);

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

	fglo = _mm_packs_epi32(fg0, fg1);
	fghi = _mm_packs_epi32(fg2, fg3);
	fg = _mm_packus_epi16(fglo, fghi);

	_mm_storeu_si128((__m128i*)FragColor, fg);
}
#endif

/////////////////////////////////////////////////////////////////////////////

void ScreenBlockDrawerSSE2::StepY()
{
	GradPosY.W += GradStepY.W;

	GradPosY.WorldX += GradStepY.WorldX;
	GradPosY.WorldY += GradStepY.WorldY;
	GradPosY.WorldZ += GradStepY.WorldZ;
	GradPosY.U += GradStepY.U;
	GradPosY.V += GradStepY.V;

	Dest += Pitch;
}

void ScreenBlockDrawerSSE2::StoreFull(int offset)
{
	_mm_storeu_si128((__m128i*)(Dest + offset), _mm_loadu_si128((const __m128i*)Shader.FragColor));
}

void ScreenBlockDrawerSSE2::StoreMasked(int offset, uint32_t mask)
{
	uint32_t *d = Dest + offset;
	for (int i = 0; i < 4; i++)
	{
		if (mask & (1 << 31))
			d[i] = Shader.FragColor[i];
		mask <<= 1;
	}
}

void ScreenBlockDrawerSSE2::ProcessMaskRange(uint32_t mask)
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

void ScreenBlockDrawerSSE2::ProcessBlock(uint32_t mask0, uint32_t mask1)
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

void ScreenBlockDrawerSSE2::SetUniforms(const TriDrawTriangleArgs *args)
{
	uint32_t maskvalue = args->uniforms->FixedLight() ? 0 : 0xffffffff;
	float *maskvaluef = (float*)&maskvalue;

	Shader.LightMask = *maskvaluef;
	Shader.Light = args->uniforms->Light() * 256.0f / 255.0f;
	Shader.Shade = 2.0f - (Shader.Light + 12.0f) / 128.0f;
	Shader.GlobVis = args->uniforms->GlobVis() * (1.0f / 32.0f);
	Shader.Light /= 256.0f;

	Shader.Tex.data = (const uint32_t *)args->uniforms->TexturePixels();
	Shader.Tex.width = args->uniforms->TextureWidth();
	Shader.Tex.height = args->uniforms->TextureHeight();
}

void ScreenBlockDrawerSSE2::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
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

void ScreenBlockDrawerSSE2::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawerSSE2 block;
	block.Dest = ((uint32_t *)args->dest) + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
