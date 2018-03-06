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

struct SWVec4f
{
	float x[4], y[4], z[4], w[4];
};

class Sampler
{
public:
	uint32_t width;
	uint32_t height;
	const uint32_t *data;

	FORCEINLINE uint32_t TextureNearest(const SWVec4f &input, int i) const;
};

class SWFragmentShader
{
public:
	// Uniforms
	Sampler Tex;
	float LightMask;
	float Light;
	float Shade;
	float GlobVis;
	PolyLight *Lights;
	int NumLights;
	float WorldNormal[3];
	uint32_t DynLightColor[4];

	// In variables
	float GradW[4];
	SWVec4f WorldPos;
	SWVec4f TexCoord;

	// Out variables
	uint32_t FragColor[4];

	FORCEINLINE void SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step);
	FORCEINLINE void Run();

private:
	FORCEINLINE uint32_t CalcDynamicLight(int i);
};

template<typename BlendT>
class ScreenBlockDrawer
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

	SWFragmentShader Shader;
};

/////////////////////////////////////////////////////////////////////////////

uint32_t Sampler::TextureNearest(const SWVec4f &input, int i) const
{
	uint32_t x = ((static_cast<uint32_t>(static_cast<int32_t>(input.x[i] * (1 << 24)) << 8) >> 16) * width) >> 16;
	uint32_t y = ((static_cast<uint32_t>(static_cast<int32_t>(input.y[i] * (1 << 24)) << 8) >> 16) * height) >> 16;
	return data[y + x * height];
}

/////////////////////////////////////////////////////////////////////////////

void SWFragmentShader::SetVaryings(ScreenTriangleStepVariables &pos, const ScreenTriangleStepVariables &step)
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

void SWFragmentShader::Run()
{
	for (int i = 0; i < 4; i++)
	{
		uint32_t fg = Tex.TextureNearest(TexCoord, i);

		float lightposf;
		if (!LightMask)
			lightposf = 1.0f - clamp(Shade - MIN(24.0f / 32.0f, GlobVis * GradW[i]), 0.0f, 31.0f / 32.0f);
		else
			lightposf = Light;

		uint32_t material_r = RPART(fg);
		uint32_t material_g = GPART(fg);
		uint32_t material_b = BPART(fg);

		uint32_t r = (int32_t)(material_r * lightposf);
		uint32_t g = (int32_t)(material_g * lightposf);
		uint32_t b = (int32_t)(material_b * lightposf);
		uint32_t a = APART(fg);

		if (NumLights != 0)
		{
			uint32_t dynlight = CalcDynamicLight(i);
			r = MIN(r + ((material_r * RPART(dynlight)) >> 8), (uint32_t)255);
			g = MIN(g + ((material_g * GPART(dynlight)) >> 8), (uint32_t)255);
			b = MIN(b + ((material_b * BPART(dynlight)) >> 8), (uint32_t)255);
		}

		FragColor[i] = MAKEARGB(a, r, g, b);
	}
}

uint32_t SWFragmentShader::CalcDynamicLight(int j)
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

	return MAKEARGB(lit[3],lit[0],lit[1],lit[2]);
}

/////////////////////////////////////////////////////////////////////////////

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::StepY()
{
	GradPosY.W += GradStepY.W;

	GradPosY.WorldX += GradStepY.WorldX;
	GradPosY.WorldY += GradStepY.WorldY;
	GradPosY.WorldZ += GradStepY.WorldZ;
	GradPosY.U += GradStepY.U;
	GradPosY.V += GradStepY.V;

	Dest += Pitch;
}

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::StoreFull(int offset)
{
	using namespace TriScreenDrawerModes;

	if (BlendT::Mode == (int)BlendModes::Opaque)
	{
		for (int i = 0; i < 4; i++)
			Dest[offset + i] = Shader.FragColor[i];
	}
	else if (BlendT::Mode == (int)BlendModes::Masked)
	{
		for (int i = 0; i < 4; i++)
		{
			if (APART(Shader.FragColor[i]))
				Dest[offset + i] = Shader.FragColor[i];
		}
	}
}

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::StoreMasked(int offset, uint32_t mask)
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

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::ProcessMaskRange(uint32_t mask)
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

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::ProcessBlock(uint32_t mask0, uint32_t mask1)
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

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::SetUniforms(const TriDrawTriangleArgs *args)
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

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::SetGradients(int destX, int destY, const ShadedTriVertex &v1, const ScreenTriangleStepVariables &gradientX, const ScreenTriangleStepVariables &gradientY)
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

template<typename BlendT>
void ScreenBlockDrawer<BlendT>::Draw(int destX, int destY, uint32_t mask0, uint32_t mask1, const TriDrawTriangleArgs *args)
{
	ScreenBlockDrawer block;
	block.Dest = ((uint32_t *)args->dest) + destX + destY * args->pitch;
	block.Pitch = args->pitch;
	block.SetUniforms(args);
	block.SetGradients(destX, destY, *args->v1, args->gradientX, args->gradientY);
	block.ProcessBlock(mask0, mask1);
}
