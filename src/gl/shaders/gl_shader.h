// 
//---------------------------------------------------------------------------
//
// Copyright(C) 2004-2016 Christoph Oelckers
// All rights reserved.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see http://www.gnu.org/licenses/
//
//--------------------------------------------------------------------------
//

#ifndef __GL_SHADERS_H__
#define __GL_SHADERS_H__

#include "gl/renderer/gl_renderstate.h"
#include "name.h"

extern bool gl_shaderactive;

enum
{
	VATTR_VERTEX = 0,
	VATTR_TEXCOORD = 1,
	VATTR_COLOR = 2,
	VATTR_VERTEX2 = 3,
	VATTR_NORMAL = 4
};

class FShaderCollection;

//==========================================================================
//
//
//==========================================================================

class FUniform1i
{
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
	}

	void Set(int newvalue)
	{
		glUniform1i(mIndex, newvalue);
	}
};

class FBufferedUniform1i
{
	int mBuffer;
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		mBuffer = 0;
	}

	void Set(int newvalue)
	{
		if (newvalue != mBuffer)
		{
			mBuffer = newvalue;
			glUniform1i(mIndex, newvalue);
		}
	}
};

class FBufferedUniform4i
{
	int mBuffer[4];
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		memset(mBuffer, 0, sizeof(mBuffer));
	}

	void Set(const int *newvalue)
	{
		if (memcmp(newvalue, mBuffer, sizeof(mBuffer)))
		{
			memcpy(mBuffer, newvalue, sizeof(mBuffer));
			glUniform4iv(mIndex, 1, newvalue);
		}
	}
};

class FBufferedUniform1f
{
	float mBuffer;
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		mBuffer = 0;
	}

	void Set(float newvalue)
	{
		if (newvalue != mBuffer)
		{
			mBuffer = newvalue;
			glUniform1f(mIndex, newvalue);
		}
	}
};

class FBufferedUniform2f
{
	float mBuffer[2];
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		memset(mBuffer, 0, sizeof(mBuffer));
	}

	void Set(const float *newvalue)
	{
		if (memcmp(newvalue, mBuffer, sizeof(mBuffer)))
		{
			memcpy(mBuffer, newvalue, sizeof(mBuffer));
			glUniform2fv(mIndex, 1, newvalue);
		}
	}

	void Set(float f1, float f2)
	{
		if (mBuffer[0] != f1 || mBuffer[1] != f2)
		{
			mBuffer[0] = f1;
			mBuffer[1] = f2;
			glUniform2fv(mIndex, 1, mBuffer);
		}
	}

};

class FBufferedUniform4f
{
	float mBuffer[4];
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		memset(mBuffer, 0, sizeof(mBuffer));
	}

	void Set(const float *newvalue)
	{
		if (memcmp(newvalue, mBuffer, sizeof(mBuffer)))
		{
			memcpy(mBuffer, newvalue, sizeof(mBuffer));
			glUniform4fv(mIndex, 1, newvalue);
		}
	}
};

class FUniform4f
{
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
	}

	void Set(const float *newvalue)
	{
		glUniform4fv(mIndex, 1, newvalue);
	}

	void Set(float a, float b, float c, float d)
	{
		glUniform4f(mIndex, a, b, c, d);
	}
};

class FBufferedUniformPE
{
	PalEntry mBuffer;
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		mBuffer = 0;
	}

	void Set(PalEntry newvalue)
	{
		if (newvalue != mBuffer)
		{
			mBuffer = newvalue;
			glUniform4f(mIndex, newvalue.r/255.f, newvalue.g/255.f, newvalue.b/255.f, newvalue.a/255.f);
		}
	}
};

class FBufferedUniformSampler
{
	int mBuffer;
	int mIndex;

public:
	void Init(GLuint hShader, const GLchar *name)
	{
		mIndex = glGetUniformLocation(hShader, name);
		mBuffer = -1;
	}

	void Set(int newvalue)
	{
		if (newvalue != mBuffer)
		{
			mBuffer = newvalue;
			glUniform1i(mIndex, newvalue);
		}
	}
};

template<typename T>
class FUsageUniform
{
public:
	void Init(GLuint hShader, const GLchar *name)
	{
		buffer.Init(hShader, name);
	}

	template<typename P1>
	void Set(P1 p1) { lastbuffer = buffer; buffer.Set(p1); TrackChanges(); }

	template<typename P1, typename P2>
	void Set(P1 p1, P2 p2) { lastbuffer = buffer; buffer.Set(p1 ,p2); TrackChanges(); }

	template<typename P1, typename P2, typename P3>
	void Set(P1 p1, P2 p2, P3 p3) { lastbuffer = buffer; buffer.Set(p1, p2, p3); TrackChanges(); }

	template<typename P1, typename P2, typename P3, typename P4>
	void Set(P1 p1, P2 p2, P3 p3, P4 p4) { lastbuffer = buffer; buffer.Set(p1, p2, p3, p4); TrackChanges(); }

	int ModifiedCount = 0;

private:
	void TrackChanges()
	{
		if (memcmp(&buffer, &lastbuffer, sizeof(T)) != 0)
			ModifiedCount++;
	}

	T buffer;
	T lastbuffer;
};


class FShader
{
	friend class FShaderCollection;
	friend class FRenderState;

	unsigned int hShader;
	unsigned int hVertProg;
	unsigned int hFragProg;
	FName mName;

	FUsageUniform<FBufferedUniform4f> muCameraPos;
	FUsageUniform<FBufferedUniform1i> muTextureMode;
	FUsageUniform<FBufferedUniform1f> muClipHeight;
	FUsageUniform<FBufferedUniform1f> muClipHeightDirection;
	FUsageUniform<FBufferedUniform2f> muClipSplit;
	FUsageUniform<FUniform4f> muClipLine;
	FUsageUniform<FBufferedUniform1f> muAlphaThreshold;

	// Colors
	FUsageUniform<FBufferedUniformPE> muObjectColor;
	FUsageUniform<FBufferedUniformPE> muObjectColor2;
	FUsageUniform<FBufferedUniform4f> muDynLightColor;
	FUsageUniform<FBufferedUniformPE> muFogColor;
	FUsageUniform<FBufferedUniform1f> muDesaturation;
	FUsageUniform<FBufferedUniform1f> muInterpolationFactor;

	// Fixed colormap stuff
	FUsageUniform<FUniform1i> muFixedColormap;
	FUsageUniform<FUniform4f> muColormapStart;
	FUsageUniform<FUniform4f> muColormapRange;

	// Glowing walls stuff
	FUsageUniform<FUniform4f> muGlowTopPlane;
	FUsageUniform<FUniform4f> muGlowTopColor;
	FUsageUniform<FUniform4f> muGlowBottomPlane;
	FUsageUniform<FUniform4f> muGlowBottomColor;

	FUsageUniform<FUniform4f> muSplitTopPlane;
	FUsageUniform<FUniform4f> muSplitBottomPlane;

	// Lighting + Fog
	FUsageUniform<FBufferedUniform4f> muLightParms;
	FUsageUniform<FBufferedUniform1i> muFogEnabled;
	FUsageUniform<FBufferedUniform1i> muPalLightLevels;
	FUsageUniform<FBufferedUniform1f> muGlobVis;

	// dynamic lights
	FUsageUniform<FBufferedUniform1i> muLightIndex;

	// Software fuzz scaling
	FUsageUniform<FBufferedUniform1i> muViewHeight;

	// Blinn glossiness and specular level
	FUsageUniform<FBufferedUniform2f> muSpecularMaterial;

	// quad drawer stuff
	// muQuadVertices
	// muQuadTexCoords
	// muQuadMode

	// matrices
	// mProjectionMatrix
	// mViewMatrix
	// mModelMatrix
	// mNormalViewMatrix
	// mNormalModelMatrix
	// mTextureMatrix

	// Timer data
	FUsageUniform<FBufferedUniform1f> muTimer;
	
	int lights_index;
	int projectionmatrix_index;
	int viewmatrix_index;
	int normalviewmatrix_index;
	int modelmatrix_index;
	int normalmodelmatrix_index;
	int texturematrix_index;
public:
	int vertexmatrix_index;
	int texcoordmatrix_index;
	int quadmode_index;
	int fakevb_index;
private:
	int currentglowstate = 0;
	int currentsplitstate = 0;
	int currentcliplinestate = 0;
	int currentfixedcolormap = 0;
	bool currentTextureMatrixState = true;// by setting the matrix state to 'true' it is guaranteed to be set the first time the render state gets applied.
	bool currentModelMatrixState = true;

public:
	FShader(const char *name)
		: mName(name)
	{
		hShader = hVertProg = hFragProg = 0;
	}

	~FShader();

	bool Load(const char * name, const char * vert_prog_lump, const char * fragprog, const char * fragprog2, const char *defines);

	void SetColormapColor(float r, float g, float b, float r1, float g1, float b1);
	void SetGlowParams(float *topcolors, float topheight, float *bottomcolors, float bottomheight);
	void SetLightRange(int start, int end, int forceadd);

	bool Bind();
	unsigned int GetHandle() const { return hShader; }

	void ApplyMatrices(VSMatrix *proj, VSMatrix *view, VSMatrix *norm);

	static FString GetUsageStats();
};

//==========================================================================
//
// The global shader manager
//
//==========================================================================
class FShaderManager
{
public:
	FShaderManager();
	~FShaderManager();

	void SetActiveShader(FShader *sh);
	FShader *GetActiveShader() const { return mActiveShader; }

	FShader *BindEffect(int effect, EPassType passType);
	FShader *Get(unsigned int eff, bool alphateston, EPassType passType);
	void ApplyMatrices(VSMatrix *proj, VSMatrix *view, EPassType passType);

	void ResetFixedColormap();

private:
	FShader *mActiveShader = nullptr;
	TArray<FShaderCollection*> mPassShaders;
};

class FShaderCollection
{
	TArray<FShader*> mMaterialShaders;
	TArray<FShader*> mMaterialShadersNAT;
	FShader *mEffectShaders[MAX_EFFECTS];

	void Clean();
	void CompileShaders(EPassType passType);
	
public:
	FShaderCollection(EPassType passType);
	~FShaderCollection();
	FShader *Compile(const char *ShaderName, const char *ShaderPath, const char *shaderdefines, bool usediscard, EPassType passType);
	int Find(const char *mame);
	FShader *BindEffect(int effect);
	void ApplyMatrices(VSMatrix *proj, VSMatrix *view);

	void ResetFixedColormap()
	{
		for (unsigned i = 0; i < mMaterialShaders.Size(); i++)
		{
			mMaterialShaders[i]->currentfixedcolormap = -1;
		}
		for (unsigned i = 0; i < mMaterialShadersNAT.Size(); i++)
		{
			mMaterialShadersNAT[i]->currentfixedcolormap = -1;
		}
	}

	FShader *Get(unsigned int eff, bool alphateston)
	{
		// indices 0-2 match the warping modes, 3 is brightmap, 4 no texture, the following are custom
		if (!alphateston && eff <= 3)
		{
			return mMaterialShadersNAT[eff];	// Non-alphatest shaders are only created for default, warp1+2 and brightmap. The rest won't get used anyway
		}
		if (eff < mMaterialShaders.Size())
		{
			return mMaterialShaders[eff];
		}
		return NULL;
	}
};

enum MaterialShaderIndex
{
	SHADER_Default,
	SHADER_Warp1,
	SHADER_Warp2,
	SHADER_Brightmap,
	SHADER_Specular,
	SHADER_SpecularBrightmap,
	SHADER_PBR,
	SHADER_PBRBrightmap,
	SHADER_NoTexture,
	SHADER_BasicFuzz,
	SHADER_SmoothFuzz,
	SHADER_SwirlyFuzz,
	SHADER_TranslucentFuzz,
	SHADER_JaggedFuzz,
	SHADER_NoiseFuzz,
	SHADER_SmoothNoiseFuzz,
	SHADER_SoftwareFuzz,
	FIRST_USER_SHADER
};

enum
{
	LIGHTBUF_BINDINGPOINT = 1
};

#endif

