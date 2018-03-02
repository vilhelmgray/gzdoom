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

class SamplerAVX2
{
public:
	uint32_t width;
	uint32_t height;
	const uint32_t *data;

	FORCEINLINE __m256i TextureNearest(const SWVec8fAVX2 &input) const;
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

	// In variables
	__m256 GradW;
	SWVec8fAVX2 WorldPos;
	SWVec8fAVX2 TexCoord;

	// Out variables
	uint32_t FragColor[8];

	FORCEINLINE void SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void Run();
};

class ScreenBlockDrawerAVX2
{
public:
	static void Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args);

private:
	FORCEINLINE void SetUniforms(const TriDrawTriangleArgs *args);
	FORCEINLINE void SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY);
	FORCEINLINE void ProcessBlock(uint32_t mask0, uint32_t mask1);
	FORCEINLINE void ProcessMaskRange(uint32_t mask);
	FORCEINLINE void StepY();
	FORCEINLINE void StoreFull();
	FORCEINLINE void StoreMasked(uint32_t mask);

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

__m256i SamplerAVX2::TextureNearest(const SWVec8fAVX2 &input) const
{
	__m256i tmpx = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.x, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmpy = _mm256_srli_epi32(_mm256_slli_epi32(_mm256_cvtps_epi32(_mm256_mul_ps(input.y, _mm256_set1_ps(1 << 24))), 8), 17);
	__m256i tmp = _mm256_packs_epi32(tmpx, tmpy);
	__m256i size = _mm256_setr_epi16(
		width << 1, width << 1, width << 1, width << 1,
		height << 1, height << 1, height << 1, height << 1,
		width << 1, width << 1, width << 1, width << 1,
		height << 1, height << 1, height << 1, height << 1);

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

	fglo = _mm256_packs_epi32(fg0, fg1);
	fghi = _mm256_packs_epi32(fg2, fg3);
	fg = _mm256_packus_epi16(fglo, fghi);

	_mm256_storeu_si256((__m256i*)FragColor, fg);
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

	Shader.Tex.data = (const uint32_t *)args->uniforms->TexturePixels();
	Shader.Tex.width = args->uniforms->TextureWidth();
	Shader.Tex.height = args->uniforms->TextureHeight();
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
