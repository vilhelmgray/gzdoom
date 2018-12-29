// 
//---------------------------------------------------------------------------
//
// Copyright(C) 2018 Magnus Norddahl
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

#ifdef WIN32
#include <d3d11.h>
#endif

#include "gl_load/gl_system.h"
#include "gl_load/gl_interface.h"
#include "gl/textures/gl_d3d11texture.h"

namespace OpenGLRenderer
{
#ifdef WIN32

	#define WGL_ACCESS_READ_ONLY_NV 0x0000
	#define WGL_ACCESS_READ_WRITE_NV 0x0001
	#define WGL_ACCESS_WRITE_DISCARD_NV 0x0002

	class GL_D3D11_TextureImpl : public GL_D3D11_Texture
	{
	public:
		typedef BOOL (*wglDXSetResourceShareHandleNVFunc)(void *dxObject, HANDLE shareHandle);
		typedef HANDLE (*wglDXOpenDeviceNVFunc)(void *dxDevice);
		typedef BOOL (*wglDXCloseDeviceNVFunc)(HANDLE hDevice);
		typedef HANDLE (*wglDXRegisterObjectNVFunc)(HANDLE hDevice, void *dxObject, GLuint name, GLenum type, GLenum access);
		typedef BOOL (*wglDXUnregisterObjectNVFunc)(HANDLE hDevice, HANDLE hObject);
		typedef BOOL (*wglDXObjectAccessNVFunc)(HANDLE hObject, GLenum access);
		typedef BOOL (*wglDXLockObjectsNVFunc)(HANDLE hDevice, GLint count, HANDLE *hObjects);
		typedef BOOL (*wglDXUnlockObjectsNVFunc)(HANDLE hDevice, GLint count, HANDLE *hObjects);

		typedef HRESULT (WINAPI *D3D11CreateDeviceFunc)(IDXGIAdapter* pAdapter, D3D_DRIVER_TYPE DriverType, HMODULE Software, UINT Flags, CONST D3D_FEATURE_LEVEL* pFeatureLevels, UINT FeatureLevels, UINT SDKVersion, ID3D11Device** ppDevice, D3D_FEATURE_LEVEL* pFeatureLevel, ID3D11DeviceContext** ppImmediateContext);

		wglDXSetResourceShareHandleNVFunc wglDXSetResourceShareHandleNV = nullptr;
		wglDXOpenDeviceNVFunc wglDXOpenDeviceNV = nullptr;
		wglDXCloseDeviceNVFunc wglDXCloseDeviceNV = nullptr;
		wglDXRegisterObjectNVFunc wglDXRegisterObjectNV = nullptr;
		wglDXUnregisterObjectNVFunc wglDXUnregisterObjectNV = nullptr;
		wglDXObjectAccessNVFunc wglDXObjectAccessNV = nullptr;
		wglDXLockObjectsNVFunc wglDXLockObjectsNV = nullptr;
		wglDXUnlockObjectsNVFunc wglDXUnlockObjectsNV = nullptr;

		HMODULE d3d11 = 0;
		D3D11CreateDeviceFunc D3D11CreateDevice = nullptr;

		ID3D11Device *device = nullptr;
		ID3D11DeviceContext *context = nullptr;
		ID3D11Texture2D *texture = nullptr;
		ID3D11Texture2D *staging = nullptr;
		HANDLE handleDevice = 0;
		HANDLE handleTexture = 0;
		D3D11_MAPPED_SUBRESOURCE mapData = { 0 };

		bool Create(GLuint gltexture, int width, int height, bool bgra)
		{
			wglDXOpenDeviceNV = (wglDXOpenDeviceNVFunc)wglGetProcAddress("wglDXOpenDeviceNV");
			if (!wglDXOpenDeviceNV)
				return false;

			wglDXCloseDeviceNV = (wglDXCloseDeviceNVFunc)wglGetProcAddress("wglDXCloseDeviceNV");
			wglDXRegisterObjectNV = (wglDXRegisterObjectNVFunc)wglGetProcAddress("wglDXRegisterObjectNV");
			wglDXUnregisterObjectNV = (wglDXUnregisterObjectNVFunc)wglGetProcAddress("wglDXUnregisterObjectNV");
			wglDXObjectAccessNV = (wglDXObjectAccessNVFunc)wglGetProcAddress("wglDXObjectAccessNV");
			wglDXLockObjectsNV = (wglDXLockObjectsNVFunc)wglGetProcAddress("wglDXLockObjectsNV");
			wglDXUnlockObjectsNV = (wglDXUnlockObjectsNVFunc)wglGetProcAddress("wglDXUnlockObjectsNV");
			wglDXSetResourceShareHandleNV = (wglDXSetResourceShareHandleNVFunc)wglGetProcAddress("wglDXSetResourceShareHandleNV");

			d3d11 = LoadLibrary("D3D11.DLL");
			if (!d3d11)
				return false;

			D3D11CreateDevice = (D3D11CreateDeviceFunc)GetProcAddress(d3d11, "D3D11CreateDevice");
			if (!D3D11CreateDevice)
				return false;

			HRESULT result = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, 0, nullptr, 0, D3D11_SDK_VERSION, &device, nullptr, &context);
			if (FAILED(result))
				return false;

			D3D11_TEXTURE2D_DESC desc = { 0 };
			desc.Width = width;
			desc.Height = height;
			desc.MipLevels = 1;
			desc.ArraySize = 1;
			desc.Format = bgra ? DXGI_FORMAT_B8G8R8A8_UNORM : DXGI_FORMAT_R8_UNORM;
			desc.SampleDesc.Count = 1;
			desc.SampleDesc.Quality = 0;
			desc.Usage = D3D11_USAGE_DEFAULT;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
			desc.CPUAccessFlags = 0;
			desc.MiscFlags = 0;
			result = device->CreateTexture2D(&desc, nullptr, &texture);
			if (FAILED(result))
				return false;

			desc.Usage = D3D11_USAGE_DYNAMIC;
			desc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
			desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
			result = device->CreateTexture2D(&desc, nullptr, &staging);
			if (FAILED(result))
				return false;

			handleDevice = wglDXOpenDeviceNV(device);
			if (handleDevice == 0)
				return false;

			handleTexture = wglDXRegisterObjectNV(handleDevice, texture, gltexture, GL_TEXTURE_2D, WGL_ACCESS_READ_ONLY_NV);
			wglDXLockObjectsNV(handleDevice, 1, &handleTexture);

			return true;
		}

		void *Map() override
		{
			wglDXUnlockObjectsNV(handleDevice, 1, &handleTexture);

			HRESULT result = context->Map(staging, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapData);
			if (FAILED(result))
			{
				mapData = { 0 };
				wglDXLockObjectsNV(handleDevice, 1, &handleTexture);
				return nullptr;
			}
			return mapData.pData;
		}

		void Unmap() override
		{
			context->Unmap(staging, 0);
			mapData = { 0 };

			context->CopyResource(texture, staging);

			context->Flush();
			wglDXLockObjectsNV(handleDevice, 1, &handleTexture);
		}

		~GL_D3D11_TextureImpl()
		{
			Cleanup();
		}

		void Cleanup()
		{
			if (handleTexture != 0)
			{
				wglDXUnlockObjectsNV(handleDevice, 1, &handleTexture);
				wglDXUnregisterObjectNV(handleDevice, handleTexture);
			}

			if (handleDevice != 0)
				wglDXCloseDeviceNV(handleDevice);

			if (staging)
				staging->Release();
			if (texture)
				texture->Release();
			if (context)
				context->Release();
			if (device)
				device->Release();

			if (d3d11)
				FreeLibrary(d3d11);

			handleTexture = 0;
			handleDevice = 0;
			staging = nullptr;
			texture = nullptr;
			context = nullptr;
			device = nullptr;
			d3d11 = 0;
		}
	};

	std::unique_ptr<GL_D3D11_Texture> GL_D3D11_Texture::Create(unsigned int *gltexture, int width, int height, bool bgra)
	{
		std::unique_ptr<GL_D3D11_TextureImpl> impl(new GL_D3D11_TextureImpl());
		glGenTextures(1, gltexture);
		if (impl->Create(*gltexture, width, height, bgra))
		{
			return impl;
		}
		else
		{
			glDeleteTextures(1, gltexture);
			return nullptr;
		}
	}

#else

	std::unique_ptr<GL_D3D11_Texture> GL_D3D11_Texture::Create(unsigned int *gltexture, int width, int height, bool bgra)
	{
		return nullptr;
	}

#endif
}
