
#pragma once

#include "gl_load/gl_interface.h"
#include <memory>

namespace OpenGLRenderer
{
	class GL_D3D11_Texture
	{
	public:
		static std::unique_ptr<GL_D3D11_Texture> Create(unsigned int *gltexture, int width, int height, bool bgra);

		virtual ~GL_D3D11_Texture() { }

		virtual void *Map() = 0;
		virtual void Unmap() = 0;
	};
}
