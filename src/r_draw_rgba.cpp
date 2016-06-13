// Emacs style mode select	 -*- C++ -*- 
//-----------------------------------------------------------------------------
//
// $Id:$
//
// Copyright (C) 1993-1996 by id Software, Inc.
//
// This source is available for distribution and/or modification
// only under the terms of the DOOM Source Code License as
// published by id Software. All rights reserved.
//
// The source is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// FITNESS FOR A PARTICULAR PURPOSE. See the DOOM Source Code License
// for more details.
//
// $Log:$
//
// DESCRIPTION:
//		True color span/column drawing functions.
//
//-----------------------------------------------------------------------------

#include <stddef.h>

#include "templates.h"
#include "doomdef.h"
#include "i_system.h"
#include "w_wad.h"
#include "r_local.h"
#include "v_video.h"
#include "doomstat.h"
#include "st_stuff.h"
#include "g_game.h"
#include "g_level.h"
#include "r_data/r_translate.h"
#include "v_palette.h"
#include "r_data/colormaps.h"
#include "r_plane.h"
#include "r_draw_rgba.h"

#include "gi.h"
#include "stats.h"
#include "x86.h"
#ifndef NO_SSE
#include <emmintrin.h>
#endif
#include <vector>

extern int vlinebits;
extern int mvlinebits;
extern int tmvlinebits;

extern "C" short spanend[MAXHEIGHT];
extern float rw_light;
extern float rw_lightstep;
extern int wallshade;

/////////////////////////////////////////////////////////////////////////////

DrawerCommandQueue *DrawerCommandQueue::Instance()
{
	static DrawerCommandQueue queue;
	return &queue;
}

DrawerCommandQueue::~DrawerCommandQueue()
{
	StopThreads();
}

void* DrawerCommandQueue::AllocMemory(size_t size)
{
	// Make sure allocations remain 16-byte aligned
	size = (size + 15) / 16 * 16;

	auto queue = Instance();
	if (queue->memorypool_pos + size > memorypool_size)
		return nullptr;

	void *data = queue->memorypool + queue->memorypool_pos;
	queue->memorypool_pos += size;
	return data;
}

void DrawerCommandQueue::Begin()
{
	auto queue = Instance();
	queue->Finish();
	queue->threaded_render++;
}

void DrawerCommandQueue::End()
{
	auto queue = Instance();
	queue->Finish();
	if (queue->threaded_render > 0)
		queue->threaded_render--;
}

void DrawerCommandQueue::WaitForWorkers()
{
	Instance()->Finish();
}

void DrawerCommandQueue::Finish()
{
	auto queue = Instance();
	if (queue->commands.empty())
		return;

	// Give worker threads something to do:

	std::unique_lock<std::mutex> start_lock(queue->start_mutex);
	queue->active_commands.swap(queue->commands);
	queue->run_id++;
	start_lock.unlock();

	queue->StartThreads();
	queue->start_condition.notify_all();

	// Do one thread ourselves:

	DrawerThread thread;
	thread.core = 0;
	thread.num_cores = queue->threads.size() + 1;

	for (int pass = 0; pass < queue->num_passes; pass++)
	{
		thread.pass_start_y = pass * queue->rows_in_pass;
		thread.pass_end_y = (pass + 1) * queue->rows_in_pass;
		if (pass + 1 == queue->num_passes)
			thread.pass_end_y = MAX(thread.pass_end_y, MAXHEIGHT);

		size_t size = queue->active_commands.size();
		for (size_t i = 0; i < size; i++)
		{
			auto &command = queue->active_commands[i];
			command->Execute(&thread);
		}
	}

	// Wait for everyone to finish:

	std::unique_lock<std::mutex> end_lock(queue->end_mutex);
	queue->end_condition.wait(end_lock, [&]() { return queue->finished_threads == queue->threads.size(); });

	// Clean up batch:

	for (auto &command : queue->active_commands)
		command->~DrawerCommand();
	queue->active_commands.clear();
	queue->memorypool_pos = 0;
	queue->finished_threads = 0;
}

void DrawerCommandQueue::StartThreads()
{
	if (!threads.empty())
		return;

	int num_threads = std::thread::hardware_concurrency();
	if (num_threads == 0)
		num_threads = 4;

	threads.resize(num_threads - 1);

	for (int i = 0; i < num_threads - 1; i++)
	{
		DrawerCommandQueue *queue = this;
		DrawerThread *thread = &threads[i];
		thread->core = i + 1;
		thread->num_cores = num_threads;
		thread->thread = std::thread([=]()
		{
			int run_id = 0;
			while (true)
			{
				// Wait until we are signalled to run:
				std::unique_lock<std::mutex> start_lock(queue->start_mutex);
				queue->start_condition.wait(start_lock, [&]() { return queue->run_id != run_id || queue->shutdown_flag; });
				if (queue->shutdown_flag)
					break;
				run_id = queue->run_id;
				start_lock.unlock();

				// Do the work:
				for (int pass = 0; pass < queue->num_passes; pass++)
				{
					thread->pass_start_y = pass * queue->rows_in_pass;
					thread->pass_end_y = (pass + 1) * queue->rows_in_pass;
					if (pass + 1 == queue->num_passes)
						thread->pass_end_y = MAX(thread->pass_end_y, MAXHEIGHT);

					size_t size = queue->active_commands.size();
					for (size_t i = 0; i < size; i++)
					{
						auto &command = queue->active_commands[i];
						command->Execute(thread);
					}
				}

				// Notify main thread that we finished:
				std::unique_lock<std::mutex> end_lock(queue->end_mutex);
				queue->finished_threads++;
				end_lock.unlock();
				queue->end_condition.notify_all();
			}
		});
	}
}

void DrawerCommandQueue::StopThreads()
{
	std::unique_lock<std::mutex> lock(start_mutex);
	shutdown_flag = true;
	lock.unlock();
	start_condition.notify_all();
	for (auto &thread : threads)
		thread.thread.join();
	threads.clear();
	lock.lock();
	shutdown_flag = false;
}

/////////////////////////////////////////////////////////////////////////////

class DrawColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_texturefrac;
	DWORD dc_iscale;
	fixed_t dc_light;
	const BYTE *dc_source;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	BYTE *dc_colormap;

public:
	DrawColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_texturefrac = ::dc_texturefrac;
		dc_iscale = ::dc_iscale;
		dc_light = ::dc_light;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		dc_colormap = ::dc_colormap;
	}

	void Execute(DrawerThread *thread) override
	{
		int 				count;
		uint32_t*			dest;
		fixed_t 			frac;
		fixed_t 			fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);

		// Zero length, column does not exceed a pixel.
		if (count <= 0)
			return;

		// Framebuffer destination address.
		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		// Determine scaling,
		//	which is the only mapping to be done.
		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		// [RH] Get local copies of these variables so that the compiler
		//		has a better chance of optimizing this well.
		const BYTE *source = dc_source;
		int pitch = dc_pitch * thread->num_cores;
		BYTE *colormap = dc_colormap;

		do
		{
			*dest = shade_pal_index(colormap[source[frac >> FRACBITS]], light, shade_constants);

			dest += pitch;
			frac += fracstep;

		} while (--count);
	}
};

class FillColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	fixed_t dc_light;
	int dc_pitch;
	int dc_color;

public:
	FillColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_light = ::dc_light;
		dc_pitch = ::dc_pitch;
		dc_color = ::dc_color;
	}

	void Execute(DrawerThread *thread) override
	{
		int 				count;
		uint32_t*			dest;

		count = thread->count_for_thread(dc_dest_y, dc_count);

		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		uint32_t light = calc_light_multiplier(dc_light);

		{
			int pitch = dc_pitch * thread->num_cores;
			uint32_t color = shade_pal_index_simple(dc_color, light);

			do
			{
				*dest = color;
				dest += pitch;
			} while (--count);
		}
	}
};

class FillAddColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	int dc_pitch;
	fixed_t dc_light;
	int dc_color;

public:
	FillAddColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_color = ::dc_color;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;

		uint32_t fg = shade_pal_index_simple(dc_color, calc_light_multiplier(dc_light));
		uint32_t fg_red = (fg >> 24) & 0xff;
		uint32_t fg_green = (fg >> 16) & 0xff;
		uint32_t fg_blue = fg & 0xff;

		do
		{
			uint32_t bg_red = (*dest >> 16) & 0xff;
			uint32_t bg_green = (*dest >> 8) & 0xff;
			uint32_t bg_blue = (*dest) & 0xff;

			uint32_t red = (fg_red + bg_red + 1) / 2;
			uint32_t green = (fg_green + bg_green + 1) / 2;
			uint32_t blue = (fg_blue + bg_blue + 1) / 2;

			*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			dest += pitch;
		} while (--count);
	}
};

class FillAddClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	int dc_pitch;
	fixed_t dc_light;
	int dc_color;

public:
	FillAddClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_color = ::dc_color;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;

		uint32_t fg = shade_pal_index_simple(dc_color, calc_light_multiplier(dc_light));
		uint32_t fg_red = (fg >> 24) & 0xff;
		uint32_t fg_green = (fg >> 16) & 0xff;
		uint32_t fg_blue = fg & 0xff;

		do
		{
			uint32_t bg_red = (*dest >> 16) & 0xff;
			uint32_t bg_green = (*dest >> 8) & 0xff;
			uint32_t bg_blue = (*dest) & 0xff;

			uint32_t red = clamp<uint32_t>(fg_red + bg_red, 0, 255);
			uint32_t green = clamp<uint32_t>(fg_green + bg_green, 0, 255);
			uint32_t blue = clamp<uint32_t>(fg_blue + bg_blue, 0, 255);

			*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			dest += pitch;
		} while (--count);
	}
};

class FillSubClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	int dc_pitch;
	int dc_color;
	fixed_t dc_light;

public:
	FillSubClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_pitch = ::dc_pitch;
		dc_color = ::dc_color;
		dc_light = ::dc_light;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;

		uint32_t fg = shade_pal_index_simple(dc_color, calc_light_multiplier(dc_light));
		uint32_t fg_red = (fg >> 24) & 0xff;
		uint32_t fg_green = (fg >> 16) & 0xff;
		uint32_t fg_blue = fg & 0xff;

		do
		{
			uint32_t bg_red = (*dest >> 16) & 0xff;
			uint32_t bg_green = (*dest >> 8) & 0xff;
			uint32_t bg_blue = (*dest) & 0xff;

			uint32_t red = clamp<uint32_t>(256 - fg_red + bg_red, 256, 256 + 255) - 255;
			uint32_t green = clamp<uint32_t>(256 - fg_green + bg_green, 256, 256 + 255) - 255;
			uint32_t blue = clamp<uint32_t>(256 - fg_blue + bg_blue, 256, 256 + 255) - 255;

			*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			dest += pitch;
		} while (--count);
	}
};

class FillRevSubClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	int dc_pitch;
	int dc_color;
	fixed_t dc_light;

public:
	FillRevSubClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_pitch = ::dc_pitch;
		dc_color = ::dc_color;
		dc_light = ::dc_light;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;

		uint32_t fg = shade_pal_index_simple(dc_color, calc_light_multiplier(dc_light));
		uint32_t fg_red = (fg >> 24) & 0xff;
		uint32_t fg_green = (fg >> 16) & 0xff;
		uint32_t fg_blue = fg & 0xff;

		do
		{
			uint32_t bg_red = (*dest >> 16) & 0xff;
			uint32_t bg_green = (*dest >> 8) & 0xff;
			uint32_t bg_blue = (*dest) & 0xff;

			uint32_t red = clamp<uint32_t>(256 + fg_red - bg_red, 256, 256 + 255) - 255;
			uint32_t green = clamp<uint32_t>(256 + fg_green - bg_green, 256, 256 + 255) - 255;
			uint32_t blue = clamp<uint32_t>(256 + fg_blue - bg_blue, 256, 256 + 255) - 255;

			*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			dest += pitch;
		} while (--count);
	}
};

class DrawFuzzColumnRGBACommand : public DrawerCommand
{
	int dc_x;
	int dc_yl;
	int dc_yh;
	BYTE *dc_destorg;
	int dc_pitch;
	int fuzzpos;
	int fuzzviewheight;

public:
	DrawFuzzColumnRGBACommand()
	{
		dc_x = ::dc_x;
		dc_yl = ::dc_yl;
		dc_yh = ::dc_yh;
		dc_destorg = ::dc_destorg;
		dc_pitch = ::dc_pitch;
		fuzzpos = ::fuzzpos;
		fuzzviewheight = ::fuzzviewheight;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;

		// Adjust borders. Low...
		if (dc_yl == 0)
			dc_yl = 1;

		// .. and high.
		if (dc_yh > fuzzviewheight)
			dc_yh = fuzzviewheight;

		count = thread->count_for_thread(dc_yl, dc_yh - dc_yl + 1);

		// Zero length.
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_yl, dc_pitch, ylookup[dc_yl] + dc_x + (uint32_t*)dc_destorg);

		int pitch = dc_pitch * thread->num_cores;
		int fuzzstep = thread->num_cores;
		int fuzz = (fuzzpos + thread->skipped_by_thread(dc_yl)) % FUZZTABLE;

		while (count > 0)
		{
			int available = (FUZZTABLE - fuzz);
			int next_wrap = available / fuzzstep;
			if (available % fuzzstep != 0)
				next_wrap++;

			int cnt = MIN(count, next_wrap);
			count -= cnt;
			do
			{
				uint32_t bg = dest[fuzzoffset[fuzz]];
				uint32_t bg_red = (bg >> 16) & 0xff;
				uint32_t bg_green = (bg >> 8) & 0xff;
				uint32_t bg_blue = (bg) & 0xff;

				uint32_t red = bg_red * 3 / 4;
				uint32_t green = bg_green * 3 / 4;
				uint32_t blue = bg_blue * 3 / 4;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				fuzz += fuzzstep;
			} while (--cnt);

			fuzz %= FUZZTABLE;
		}
	}
};

class DrawAddColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	BYTE *dc_colormap;

public:
	DrawAddColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		dc_colormap = ::dc_colormap;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;

			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;
			BYTE *colormap = dc_colormap;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(colormap[source[frac >> FRACBITS]], light, shade_constants);

				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawTranslatedColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	BYTE *dc_translation;
	const BYTE *dc_source;
	int dc_pitch;

public:
	DrawTranslatedColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_translation = ::dc_translation;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
	}

	void Execute(DrawerThread *thread) override
	{
		int 				count;
		uint32_t*			dest;
		fixed_t 			frac;
		fixed_t 			fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			// [RH] Local copies of global vars to improve compiler optimizations
			BYTE *translation = dc_translation;
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;

			do
			{
				*dest = shade_pal_index(translation[source[frac >> FRACBITS]], light, shade_constants);
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawTlatedAddColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	BYTE *dc_translation;
	const BYTE *dc_source;
	int dc_pitch;

public:
	DrawTlatedAddColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_translation = ::dc_translation;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			BYTE *translation = dc_translation;
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(translation[source[frac >> FRACBITS]], light, shade_constants);

				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawShadedColumnRGBACommand : public DrawerCommand
{
private:
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	fixed_t dc_light;
	const BYTE *dc_source;
	lighttable_t *dc_colormap;
	int dc_color;
	int dc_pitch;

public:
	DrawShadedColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_light = ::dc_light;
		dc_source = ::dc_source;
		dc_colormap = ::dc_colormap;
		dc_color = ::dc_color;
		dc_pitch = ::dc_pitch;
	}

	void Execute(DrawerThread *thread) override
	{
		int  count;
		uint32_t *dest;
		fixed_t frac, fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);

		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		uint32_t fg = shade_pal_index_simple(dc_color, calc_light_multiplier(dc_light));
		uint32_t fg_red = (fg >> 16) & 0xff;
		uint32_t fg_green = (fg >> 8) & 0xff;
		uint32_t fg_blue = fg & 0xff;

		{
			const BYTE *source = dc_source;
			BYTE *colormap = dc_colormap;
			int pitch = dc_pitch * thread->num_cores;

			do
			{
				DWORD alpha = clamp<DWORD>(colormap[source[frac >> FRACBITS]], 0, 64);
				DWORD inv_alpha = 64 - alpha;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = (fg_red * alpha + bg_red * inv_alpha) / 64;
				uint32_t green = (fg_green * alpha + bg_green * inv_alpha) / 64;
				uint32_t blue = (fg_blue * alpha + bg_blue * inv_alpha) / 64;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawAddClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	DrawAddClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(source[frac >> FRACBITS], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawAddClampTranslatedColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	BYTE *dc_translation;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	DrawAddClampTranslatedColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_translation = ::dc_translation;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			BYTE *translation = dc_translation;
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(translation[source[frac >> FRACBITS]], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawSubClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	DrawSubClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(source[frac >> FRACBITS], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 - fg_red * fg_alpha + bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 - fg_green * fg_alpha + bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 - fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawSubClampTranslatedColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	BYTE *dc_translation;

public:
	DrawSubClampTranslatedColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		dc_translation = ::dc_translation;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			BYTE *translation = dc_translation;
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(translation[source[frac >> FRACBITS]], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 - fg_red * fg_alpha + bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 - fg_green * fg_alpha + bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 - fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawRevSubClampColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	DrawRevSubClampColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;
			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(source[frac >> FRACBITS], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 + fg_red * fg_alpha - bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 + fg_green * fg_alpha - bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 + fg_blue * fg_alpha - bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawRevSubClampTranslatedColumnRGBACommand : public DrawerCommand
{
	int dc_count;
	BYTE *dc_dest;
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	const BYTE *dc_source;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	BYTE *dc_translation;

public:
	DrawRevSubClampTranslatedColumnRGBACommand()
	{
		dc_count = ::dc_count;
		dc_dest = ::dc_dest;
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_source = ::dc_source;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		dc_translation = ::dc_translation;
	}

	void Execute(DrawerThread *thread) override
	{
		int count;
		uint32_t *dest;
		fixed_t frac;
		fixed_t fracstep;

		count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);

		fracstep = dc_iscale * thread->num_cores;
		frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);

		{
			BYTE *translation = dc_translation;
			const BYTE *source = dc_source;
			int pitch = dc_pitch * thread->num_cores;
			uint32_t light = calc_light_multiplier(dc_light);
			ShadeConstants shade_constants = dc_shade_constants;

			uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
			uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

			do
			{
				uint32_t fg = shade_pal_index(translation[source[frac >> FRACBITS]], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 + fg_red * fg_alpha - bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 + fg_green * fg_alpha - bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 + fg_blue * fg_alpha - bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				dest += pitch;
				frac += fracstep;
			} while (--count);
		}
	}
};

class DrawSpanRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_x1;
	int ds_x2;
	int ds_y;
	int ds_xbits;
	int ds_ybits;
	BYTE *dc_destorg;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;

public:
	DrawSpanRGBACommand()
	{
		ds_source = (const uint32_t*)::ds_source;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
		dc_destorg = ::dc_destorg;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
	}

#ifdef NO_SSE
	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.

			do
			{
				// Current texture index in u,v.
				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));

				// Lookup pixel from flat texture tile
				*dest++ = shade_bgra(source[spot], light, shade_constants);

				// Next step in u,v.
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;

			do
			{
				// Current texture index in u,v.
				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);

				// Lookup pixel from flat texture tile
				*dest++ = shade_bgra(source[spot], light, shade_constants);

				// Next step in u,v.
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
#else
	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.

			int sse_count = count / 4;
			count -= sse_count * 4;

			if (shade_constants.simple_shade)
			{
				SSE_SHADE_SIMPLE_INIT(light);

				while (sse_count--)
				{
					// Current texture index in u,v.
					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p0 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p1 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p2 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p3 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					// Lookup pixel from flat texture tile,
					//  re-index using light/colormap.
					__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
					SSE_SHADE_SIMPLE(fg);
					_mm_storeu_si128((__m128i*)dest, fg);

					// Next step in u,v.
					dest += 4;
				}
			}
			else
			{
				SSE_SHADE_INIT(light, shade_constants);

				while (sse_count--)
				{
					// Current texture index in u,v.
					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p0 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p1 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p2 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
					uint32_t p3 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					// Lookup pixel from flat texture tile,
					//  re-index using light/colormap.
					__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
					SSE_SHADE(fg, shade_constants);
					_mm_storeu_si128((__m128i*)dest, fg);

					// Next step in u,v.
					dest += 4;
				}
			}

			if (count == 0)
				return;

			do
			{
				// Current texture index in u,v.
				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));

				// Lookup pixel from flat texture tile
				*dest++ = shade_bgra(source[spot], light, shade_constants);

				// Next step in u,v.
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;

			int sse_count = count / 4;
			count -= sse_count * 4;

			if (shade_constants.simple_shade)
			{
				SSE_SHADE_SIMPLE_INIT(light);

				while (sse_count--)
				{
					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p0 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p1 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p2 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p3 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					// Lookup pixel from flat texture tile
					__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
					SSE_SHADE_SIMPLE(fg);
					_mm_storeu_si128((__m128i*)dest, fg);
					dest += 4;
				}
			}
			else
			{
				SSE_SHADE_INIT(light, shade_constants);

				while (sse_count--)
				{
					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p0 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p1 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p2 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
					uint32_t p3 = source[spot];
					xfrac += xstep;
					yfrac += ystep;

					// Lookup pixel from flat texture tile
					__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
					SSE_SHADE(fg, shade_constants);
					_mm_storeu_si128((__m128i*)dest, fg);
					dest += 4;
				}
			}

			if (count == 0)
				return;

			do
			{
				// Current texture index in u,v.
				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);

				// Lookup pixel from flat texture tile
				*dest++ = shade_bgra(source[spot], light, shade_constants);

				// Next step in u,v.
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
#endif
};

class DrawSpanMaskedRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	BYTE *dc_destorg;
	int ds_x1;
	int ds_x2;
	int ds_y1;
	int ds_y;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_xbits;
	int ds_ybits;

public:
	DrawSpanMaskedRGBACommand()
	{
		ds_source = (const uint32_t*)::ds_source;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		dc_destorg = ::dc_destorg;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
				texdata = source[spot];
				if (texdata != 0)
				{
					*dest = shade_bgra(texdata, light, shade_constants);
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
				texdata = source[spot];
				if (texdata != 0)
				{
					*dest = shade_bgra(texdata, light, shade_constants);
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
};

class DrawSpanTranslucentRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	BYTE *dc_destorg;
	int ds_x1;
	int ds_x2;
	int ds_y1;
	int ds_y;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_xbits;
	int ds_ybits;

public:
	DrawSpanTranslucentRGBACommand()
	{
		ds_source = (const uint32_t *)::ds_source;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		dc_destorg = ::dc_destorg;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.
			do
			{
				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));

				uint32_t fg = shade_bgra(source[spot], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = (fg) & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
				uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
				uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

				*dest++ = 0xff000000 | (red << 16) | (green << 8) | blue;

				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;
			do
			{
				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);

				uint32_t fg = shade_bgra(source[spot], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = (fg) & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
				uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
				uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

				*dest++ = 0xff000000 | (red << 16) | (green << 8) | blue;

				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
};

class DrawSpanMaskedTranslucentRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	BYTE *dc_destorg;
	int ds_x1;
	int ds_x2;
	int ds_y1;
	int ds_y;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_xbits;
	int ds_ybits;

public:
	DrawSpanMaskedTranslucentRGBACommand()
	{
		ds_source = (const uint32_t*)::ds_source;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		dc_destorg = ::dc_destorg;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
				texdata = source[spot];
				if (texdata != 0)
				{
					uint32_t fg = shade_bgra(texdata, light, shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = (fg) & 0xff;

					uint32_t bg_red = (*dest >> 16) & 0xff;
					uint32_t bg_green = (*dest >> 8) & 0xff;
					uint32_t bg_blue = (*dest) & 0xff;

					uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
					uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
					uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

					*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
				texdata = source[spot];
				if (texdata != 0)
				{
					uint32_t fg = shade_bgra(texdata, light, shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = (fg) & 0xff;

					uint32_t bg_red = (*dest >> 16) & 0xff;
					uint32_t bg_green = (*dest >> 8) & 0xff;
					uint32_t bg_blue = (*dest) & 0xff;

					uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
					uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
					uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

					*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
};

class DrawSpanAddClampRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	BYTE *dc_destorg;
	int ds_x1;
	int ds_x2;
	int ds_y1;
	int ds_y;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_xbits;
	int ds_ybits;

public:
	DrawSpanAddClampRGBACommand()
	{
		ds_source = (const uint32_t*)::ds_source;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		dc_destorg = ::dc_destorg;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.
			do
			{
				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));

				uint32_t fg = shade_bgra(source[spot], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = (fg) & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest++ = 0xff000000 | (red << 16) | (green << 8) | blue;

				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;
			do
			{
				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);

				uint32_t fg = shade_bgra(source[spot], light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = (fg) & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest++ = 0xff000000 | (red << 16) | (green << 8) | blue;

				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
};

class DrawSpanMaskedAddClampRGBACommand : public DrawerCommand
{
	const uint32_t *ds_source;
	fixed_t ds_light;
	ShadeConstants ds_shade_constants;
	fixed_t ds_xfrac;
	fixed_t ds_yfrac;
	BYTE *dc_destorg;
	int ds_x1;
	int ds_x2;
	int ds_y1;
	int ds_y;
	fixed_t ds_xstep;
	fixed_t ds_ystep;
	int ds_xbits;
	int ds_ybits;

public:
	DrawSpanMaskedAddClampRGBACommand()
	{
		ds_source = (const uint32_t*)::ds_source;
		ds_light = ::ds_light;
		ds_shade_constants = ::ds_shade_constants;
		ds_xfrac = ::ds_xfrac;
		ds_yfrac = ::ds_yfrac;
		dc_destorg = ::dc_destorg;
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		ds_xstep = ::ds_xstep;
		ds_ystep = ::ds_ystep;
		ds_xbits = ::ds_xbits;
		ds_ybits = ::ds_ybits;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		dsfixed_t			xfrac;
		dsfixed_t			yfrac;
		dsfixed_t			xstep;
		dsfixed_t			ystep;
		uint32_t*			dest;
		const uint32_t*		source = ds_source;
		int 				count;
		int 				spot;

		uint32_t light = calc_light_multiplier(ds_light);
		ShadeConstants shade_constants = ds_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		xfrac = ds_xfrac;
		yfrac = ds_yfrac;

		dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;

		count = ds_x2 - ds_x1 + 1;

		xstep = ds_xstep;
		ystep = ds_ystep;

		if (ds_xbits == 6 && ds_ybits == 6)
		{
			// 64x64 is the most common case by far, so special case it.
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> (32 - 6 - 6))&(63 * 64)) + (yfrac >> (32 - 6));
				texdata = source[spot];
				if (texdata != 0)
				{
					uint32_t fg = shade_bgra(texdata, light, shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = (fg) & 0xff;

					uint32_t bg_red = (*dest >> 16) & 0xff;
					uint32_t bg_green = (*dest >> 8) & 0xff;
					uint32_t bg_blue = (*dest) & 0xff;

					uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
					uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
					uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

					*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
		else
		{
			BYTE yshift = 32 - ds_ybits;
			BYTE xshift = yshift - ds_xbits;
			int xmask = ((1 << ds_xbits) - 1) << ds_ybits;
			do
			{
				uint32_t texdata;

				spot = ((xfrac >> xshift) & xmask) + (yfrac >> yshift);
				texdata = source[spot];
				if (texdata != 0)
				{
					uint32_t fg = shade_bgra(texdata, light, shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = (fg) & 0xff;

					uint32_t bg_red = (*dest >> 16) & 0xff;
					uint32_t bg_green = (*dest >> 8) & 0xff;
					uint32_t bg_blue = (*dest) & 0xff;

					uint32_t red = (fg_red * fg_alpha + bg_red * bg_alpha) / 256;
					uint32_t green = (fg_green * fg_alpha + bg_green * bg_alpha) / 256;
					uint32_t blue = (fg_blue * fg_alpha + bg_blue * bg_alpha) / 256;

					*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				dest++;
				xfrac += xstep;
				yfrac += ystep;
			} while (--count);
		}
	}
};

class FillSpanRGBACommand : public DrawerCommand
{
	int ds_x1;
	int ds_x2;
	int ds_y;
	BYTE *dc_destorg;
	fixed_t ds_light;
	int ds_color;

public:
	FillSpanRGBACommand()
	{
		ds_x1 = ::ds_x1;
		ds_x2 = ::ds_x2;
		ds_y = ::ds_y;
		dc_destorg = ::dc_destorg;
		ds_light = ::ds_light;
		ds_color = ::ds_color;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(ds_y))
			return;

		uint32_t *dest = ylookup[ds_y] + ds_x1 + (uint32_t*)dc_destorg;
		int count = (ds_x2 - ds_x1 + 1);
		uint32_t light = calc_light_multiplier(ds_light);
		uint32_t color = shade_pal_index_simple(ds_color, light);
		for (int i = 0; i < count; i++)
			dest[i] = color;
	}
};

class Vlinec1RGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int vlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;

public:
	Vlinec1RGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		vlinebits = ::vlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = vlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		do
		{
			*dest = shade_bgra(source[frac >> bits], light, shade_constants);
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Vlinec4RGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	int vlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Vlinec4RGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		vlinebits = ::vlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

#ifdef NO_SSE
	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = vlinebits;
		DWORD place;

		uint32_t light0 = calc_light_multiplier(palookuplight[0]);
		uint32_t light1 = calc_light_multiplier(palookuplight[1]);
		uint32_t light2 = calc_light_multiplier(palookuplight[2]);
		uint32_t light3 = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			dest[0] = shade_bgra(bufplce[0][(place = local_vplce[0]) >> bits], light0, shade_constants); local_vplce[0] = place + local_vince[0];
			dest[1] = shade_bgra(bufplce[1][(place = local_vplce[1]) >> bits], light1, shade_constants); local_vplce[1] = place + local_vince[1];
			dest[2] = shade_bgra(bufplce[2][(place = local_vplce[2]) >> bits], light2, shade_constants); local_vplce[2] = place + local_vince[2];
			dest[3] = shade_bgra(bufplce[3][(place = local_vplce[3]) >> bits], light3, shade_constants); local_vplce[3] = place + local_vince[3];
			dest += pitch;
		} while (--count);
	}
#else
	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = vlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light0 = calc_light_multiplier(palookuplight[0]);
		uint32_t light1 = calc_light_multiplier(palookuplight[1]);
		uint32_t light2 = calc_light_multiplier(palookuplight[2]);
		uint32_t light3 = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		if (shade_constants.simple_shade)
		{
			SSE_SHADE_SIMPLE_INIT4(light3, light2, light1, light0);
			do
			{
				DWORD place0 = local_vplce[0];
				DWORD place1 = local_vplce[1];
				DWORD place2 = local_vplce[2];
				DWORD place3 = local_vplce[3];

				uint32_t p0 = bufplce[0][place0 >> bits];
				uint32_t p1 = bufplce[1][place1 >> bits];
				uint32_t p2 = bufplce[2][place2 >> bits];
				uint32_t p3 = bufplce[3][place3 >> bits];

				local_vplce[0] = place0 + local_vince[0];
				local_vplce[1] = place1 + local_vince[1];
				local_vplce[2] = place2 + local_vince[2];
				local_vplce[3] = place3 + local_vince[3];

				__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
				SSE_SHADE_SIMPLE(fg);
				_mm_storeu_si128((__m128i*)dest, fg);
				dest += pitch;
			} while (--count);
		}
		else
		{
			SSE_SHADE_INIT4(light3, light2, light1, light0, shade_constants);
			do
			{
				DWORD place0 = local_vplce[0];
				DWORD place1 = local_vplce[1];
				DWORD place2 = local_vplce[2];
				DWORD place3 = local_vplce[3];

				uint32_t p0 = bufplce[0][place0 >> bits];
				uint32_t p1 = bufplce[1][place1 >> bits];
				uint32_t p2 = bufplce[2][place2 >> bits];
				uint32_t p3 = bufplce[3][place3 >> bits];

				local_vplce[0] = place0 + local_vince[0];
				local_vplce[1] = place1 + local_vince[1];
				local_vplce[2] = place2 + local_vince[2];
				local_vplce[3] = place3 + local_vince[3];

				__m128i fg = _mm_set_epi32(p3, p2, p1, p0);
				SSE_SHADE(fg, shade_constants);
				_mm_storeu_si128((__m128i*)dest, fg);
				dest += pitch;
			} while (--count);
		}
	}
#endif
};

class Mvlinec1RGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int mvlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;

public:
	Mvlinec1RGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		mvlinebits = ::mvlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = mvlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		do
		{
			uint32_t pix = source[frac >> bits];
			if (pix != 0)
			{
				*dest = shade_bgra(pix, light, shade_constants);
			}
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Mvlinec4RGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	int mvlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Mvlinec4RGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		mvlinebits = ::mvlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

#ifdef NO_SSE
	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = mvlinebits;
		DWORD place;

		uint32_t light0 = calc_light_multiplier(palookuplight[0]);
		uint32_t light1 = calc_light_multiplier(palookuplight[1]);
		uint32_t light2 = calc_light_multiplier(palookuplight[2]);
		uint32_t light3 = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			uint32_t pix;
			pix = bufplce[0][(place = local_vplce[0]) >> bits]; if (pix) dest[0] = shade_bgra(pix, light0, shade_constants); local_vplce[0] = place + local_vince[0];
			pix = bufplce[1][(place = local_vplce[1]) >> bits]; if (pix) dest[1] = shade_bgra(pix, light1, shade_constants); local_vplce[1] = place + local_vince[1];
			pix = bufplce[2][(place = local_vplce[2]) >> bits]; if (pix) dest[2] = shade_bgra(pix, light2, shade_constants); local_vplce[2] = place + local_vince[2];
			pix = bufplce[3][(place = local_vplce[3]) >> bits]; if (pix) dest[3] = shade_bgra(pix, light3, shade_constants); local_vplce[3] = place + local_vince[3];
			dest += pitch;
		} while (--count);
	}
#else
	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = mvlinebits;

		uint32_t light0 = calc_light_multiplier(palookuplight[0]);
		uint32_t light1 = calc_light_multiplier(palookuplight[1]);
		uint32_t light2 = calc_light_multiplier(palookuplight[2]);
		uint32_t light3 = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		if (shade_constants.simple_shade)
		{
			SSE_SHADE_SIMPLE_INIT4(light3, light2, light1, light0);
			do
			{
				DWORD place0 = local_vplce[0];
				DWORD place1 = local_vplce[1];
				DWORD place2 = local_vplce[2];
				DWORD place3 = local_vplce[3];

				uint32_t pix0 = bufplce[0][place0 >> bits];
				uint32_t pix1 = bufplce[1][place1 >> bits];
				uint32_t pix2 = bufplce[2][place2 >> bits];
				uint32_t pix3 = bufplce[3][place3 >> bits];

				// movemask = !(pix == 0)
				__m128i movemask = _mm_xor_si128(_mm_cmpeq_epi32(_mm_set_epi32(pix3, pix2, pix1, pix0), _mm_setzero_si128()), _mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128()));

				local_vplce[0] = place0 + local_vince[0];
				local_vplce[1] = place1 + local_vince[1];
				local_vplce[2] = place2 + local_vince[2];
				local_vplce[3] = place3 + local_vince[3];

				__m128i fg = _mm_set_epi32(pix3, pix2, pix1, pix0);
				SSE_SHADE_SIMPLE(fg);
				_mm_maskmoveu_si128(fg, movemask, (char*)dest);
				dest += pitch;
			} while (--count);
		}
		else
		{
			SSE_SHADE_INIT4(light3, light2, light1, light0, shade_constants);
			do
			{
				DWORD place0 = local_vplce[0];
				DWORD place1 = local_vplce[1];
				DWORD place2 = local_vplce[2];
				DWORD place3 = local_vplce[3];

				uint32_t pix0 = bufplce[0][place0 >> bits];
				uint32_t pix1 = bufplce[1][place1 >> bits];
				uint32_t pix2 = bufplce[2][place2 >> bits];
				uint32_t pix3 = bufplce[3][place3 >> bits];

				// movemask = !(pix == 0)
				__m128i movemask = _mm_xor_si128(_mm_cmpeq_epi32(_mm_set_epi32(pix3, pix2, pix1, pix0), _mm_setzero_si128()), _mm_cmpeq_epi32(_mm_setzero_si128(), _mm_setzero_si128()));

				local_vplce[0] = place0 + local_vince[0];
				local_vplce[1] = place1 + local_vince[1];
				local_vplce[2] = place2 + local_vince[2];
				local_vplce[3] = place3 + local_vince[3];

				__m128i fg = _mm_set_epi32(pix3, pix2, pix1, pix0);
				SSE_SHADE(fg, shade_constants);
				_mm_maskmoveu_si128(fg, movemask, (char*)dest);
				dest += pitch;
			} while (--count);
		}
	}
#endif
};

class Tmvline1AddRGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int tmvlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	Tmvline1AddRGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		tmvlinebits = ::tmvlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = tmvlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		do
		{
			uint32_t pix = source[frac >> bits];
			if (pix != 0)
			{
				uint32_t fg = shade_bgra(pix, light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			}
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Tmvline4AddRGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	int tmvlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Tmvline4AddRGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		tmvlinebits = ::tmvlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = tmvlinebits;

		uint32_t light[4];
		light[0] = calc_light_multiplier(palookuplight[0]);
		light[1] = calc_light_multiplier(palookuplight[1]);
		light[2] = calc_light_multiplier(palookuplight[2]);
		light[3] = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			for (int i = 0; i < 4; ++i)
			{
				uint32_t pix = bufplce[i][local_vplce[i] >> bits];
				if (pix != 0)
				{
					uint32_t fg = shade_bgra(pix, light[i], shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = fg & 0xff;

					uint32_t bg_red = (*dest >> 16) & 0xff;
					uint32_t bg_green = (*dest >> 8) & 0xff;
					uint32_t bg_blue = (*dest) & 0xff;

					uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
					uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
					uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

					dest[i] = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				local_vplce[i] += local_vince[i];
			}
			dest += pitch;
		} while (--count);
	}
};

class Tmvline1AddClampRGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int tmvlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	Tmvline1AddClampRGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		tmvlinebits = ::tmvlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = tmvlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		do
		{
			uint32_t pix = source[frac >> bits];
			if (pix != 0)
			{
				uint32_t fg = shade_bgra(pix, light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
				uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
				uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			}
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Tmvline4AddClampRGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	int tmvlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Tmvline4AddClampRGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		tmvlinebits = ::tmvlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = tmvlinebits;

		uint32_t light[4];
		light[0] = calc_light_multiplier(palookuplight[0]);
		light[1] = calc_light_multiplier(palookuplight[1]);
		light[2] = calc_light_multiplier(palookuplight[2]);
		light[3] = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			for (int i = 0; i < 4; ++i)
			{
				uint32_t pix = bufplce[i][local_vplce[i] >> bits];
				if (pix != 0)
				{
					uint32_t fg = shade_bgra(pix, light[i], shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = fg & 0xff;

					uint32_t bg_red = (dest[i] >> 16) & 0xff;
					uint32_t bg_green = (dest[i] >> 8) & 0xff;
					uint32_t bg_blue = (dest[i]) & 0xff;

					uint32_t red = clamp<uint32_t>((fg_red * fg_alpha + bg_red * bg_alpha) / 256, 0, 255);
					uint32_t green = clamp<uint32_t>((fg_green * fg_alpha + bg_green * bg_alpha) / 256, 0, 255);
					uint32_t blue = clamp<uint32_t>((fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 0, 255);

					dest[i] = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				local_vplce[i] += local_vince[i];
			}
			dest += pitch;
		} while (--count);
	}
};

class Tmvline1SubClampRGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int tmvlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	Tmvline1SubClampRGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		tmvlinebits = ::tmvlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = tmvlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		do
		{
			uint32_t pix = source[frac >> bits];
			if (pix != 0)
			{
				uint32_t fg = shade_bgra(pix, light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 - fg_red * fg_alpha + bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 - fg_green * fg_alpha + bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 - fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			}
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Tmvline4SubClampRGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	int tmvlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Tmvline4SubClampRGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		tmvlinebits = ::tmvlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = tmvlinebits;

		uint32_t light[4];
		light[0] = calc_light_multiplier(palookuplight[0]);
		light[1] = calc_light_multiplier(palookuplight[1]);
		light[2] = calc_light_multiplier(palookuplight[2]);
		light[3] = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			for (int i = 0; i < 4; ++i)
			{
				uint32_t pix = bufplce[i][local_vplce[i] >> bits];
				if (pix != 0)
				{
					uint32_t fg = shade_bgra(pix, light[i], shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = fg & 0xff;

					uint32_t bg_red = (dest[i] >> 16) & 0xff;
					uint32_t bg_green = (dest[i] >> 8) & 0xff;
					uint32_t bg_blue = (dest[i]) & 0xff;

					uint32_t red = clamp<uint32_t>((0x10000 - fg_red * fg_alpha + bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
					uint32_t green = clamp<uint32_t>((0x10000 - fg_green * fg_alpha + bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
					uint32_t blue = clamp<uint32_t>((0x10000 - fg_blue * fg_alpha + bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

					dest[i] = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				local_vplce[i] += local_vince[i];
			}
			dest += pitch;
		} while (--count);
	}
};

class Tmvline1RevSubClampRGBACommand : public DrawerCommand
{
	DWORD dc_iscale;
	DWORD dc_texturefrac;
	int dc_count;
	const BYTE *dc_source;
	BYTE *dc_dest;
	int tmvlinebits;
	int dc_pitch;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;

public:
	Tmvline1RevSubClampRGBACommand()
	{
		dc_iscale = ::dc_iscale;
		dc_texturefrac = ::dc_texturefrac;
		dc_count = ::dc_count;
		dc_source = ::dc_source;
		dc_dest = ::dc_dest;
		tmvlinebits = ::tmvlinebits;
		dc_pitch = ::dc_pitch;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		DWORD fracstep = dc_iscale * thread->num_cores;
		DWORD frac = dc_texturefrac + dc_iscale * thread->skipped_by_thread(dc_dest_y);
		const uint32 *source = (const uint32 *)dc_source;
		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int bits = tmvlinebits;
		int pitch = dc_pitch * thread->num_cores;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		do
		{
			uint32_t pix = source[frac >> bits];
			if (pix != 0)
			{
				uint32_t fg = shade_bgra(pix, light, shade_constants);
				uint32_t fg_red = (fg >> 16) & 0xff;
				uint32_t fg_green = (fg >> 8) & 0xff;
				uint32_t fg_blue = fg & 0xff;

				uint32_t bg_red = (*dest >> 16) & 0xff;
				uint32_t bg_green = (*dest >> 8) & 0xff;
				uint32_t bg_blue = (*dest) & 0xff;

				uint32_t red = clamp<uint32_t>((0x10000 + fg_red * fg_alpha - bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t green = clamp<uint32_t>((0x10000 + fg_green * fg_alpha - bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
				uint32_t blue = clamp<uint32_t>((0x10000 + fg_blue * fg_alpha - bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

				*dest = 0xff000000 | (red << 16) | (green << 8) | blue;
			}
			frac += fracstep;
			dest += pitch;
		} while (--count);
	}
};

class Tmvline4RevSubClampRGBACommand : public DrawerCommand
{
	BYTE *dc_dest;
	int dc_count;
	int dc_pitch;
	ShadeConstants dc_shade_constants;
	fixed_t dc_srcalpha;
	fixed_t dc_destalpha;
	int tmvlinebits;
	fixed_t palookuplight[4];
	DWORD vplce[4];
	DWORD vince[4];
	const uint32 *bufplce[4];

public:
	Tmvline4RevSubClampRGBACommand()
	{
		dc_dest = ::dc_dest;
		dc_count = ::dc_count;
		dc_pitch = ::dc_pitch;
		dc_shade_constants = ::dc_shade_constants;
		dc_srcalpha = ::dc_srcalpha;
		dc_destalpha = ::dc_destalpha;
		tmvlinebits = ::tmvlinebits;
		for (int i = 0; i < 4; i++)
		{
			palookuplight[i] = ::palookuplight[i];
			vplce[i] = ::vplce[i];
			vince[i] = ::vince[i];
			bufplce[i] = (const uint32 *)::bufplce[i];
		}
	}

	void Execute(DrawerThread *thread) override
	{
		int count = thread->count_for_thread(dc_dest_y, dc_count);
		if (count <= 0)
			return;

		uint32_t *dest = thread->dest_for_thread(dc_dest_y, dc_pitch, (uint32_t*)dc_dest);
		int pitch = dc_pitch * thread->num_cores;
		int bits = tmvlinebits;

		uint32_t light[4];
		light[0] = calc_light_multiplier(palookuplight[0]);
		light[1] = calc_light_multiplier(palookuplight[1]);
		light[2] = calc_light_multiplier(palookuplight[2]);
		light[3] = calc_light_multiplier(palookuplight[3]);

		ShadeConstants shade_constants = dc_shade_constants;

		uint32_t fg_alpha = dc_srcalpha >> (FRACBITS - 8);
		uint32_t bg_alpha = dc_destalpha >> (FRACBITS - 8);

		DWORD local_vplce[4] = { vplce[0], vplce[1], vplce[2], vplce[3] };
		DWORD local_vince[4] = { vince[0], vince[1], vince[2], vince[3] };
		int skipped = thread->skipped_by_thread(dc_dest_y);
		for (int i = 0; i < 4; i++)
		{
			local_vplce[i] += local_vince[i] * skipped;
			local_vince[i] *= thread->num_cores;
		}

		do
		{
			for (int i = 0; i < 4; ++i)
			{
				uint32_t pix = bufplce[i][local_vplce[i] >> bits];
				if (pix != 0)
				{
					uint32_t fg = shade_bgra(pix, light[i], shade_constants);
					uint32_t fg_red = (fg >> 16) & 0xff;
					uint32_t fg_green = (fg >> 8) & 0xff;
					uint32_t fg_blue = fg & 0xff;

					uint32_t bg_red = (dest[i] >> 16) & 0xff;
					uint32_t bg_green = (dest[i] >> 8) & 0xff;
					uint32_t bg_blue = (dest[i]) & 0xff;

					uint32_t red = clamp<uint32_t>((0x10000 + fg_red * fg_alpha - bg_red * bg_alpha) / 256, 256, 256 + 255) - 256;
					uint32_t green = clamp<uint32_t>((0x10000 + fg_green * fg_alpha - bg_green * bg_alpha) / 256, 256, 256 + 255) - 256;
					uint32_t blue = clamp<uint32_t>((0x10000 + fg_blue * fg_alpha - bg_blue * bg_alpha) / 256, 256, 256 + 255) - 256;

					dest[i] = 0xff000000 | (red << 16) | (green << 8) | blue;
				}
				local_vplce[i] += local_vince[i];
			}
			dest += pitch;
		} while (--count);
	}
};

class DrawFogBoundaryLineRGBACommand : public DrawerCommand
{
	int _y;
	int _x;
	int _x2;
	BYTE *dc_destorg;
	fixed_t dc_light;
	ShadeConstants dc_shade_constants;

public:
	DrawFogBoundaryLineRGBACommand(int y, int x, int x2)
	{
		_y = y;
		_x = x;
		_x2 = x2;

		dc_destorg = ::dc_destorg;
		dc_light = ::dc_light;
		dc_shade_constants = ::dc_shade_constants;
	}

	void Execute(DrawerThread *thread) override
	{
		if (thread->line_skipped_by_thread(_y))
			return;

		int y = _y;
		int x = _x;
		int x2 = _x2;

		uint32_t *dest = ylookup[y] + (uint32_t*)dc_destorg;

		uint32_t light = calc_light_multiplier(dc_light);
		ShadeConstants constants = dc_shade_constants;

		do
		{
			uint32_t red = (dest[x] >> 16) & 0xff;
			uint32_t green = (dest[x] >> 8) & 0xff;
			uint32_t blue = dest[x] & 0xff;

			if (constants.simple_shade)
			{
				red = red * light / 256;
				green = green * light / 256;
				blue = blue * light / 256;
			}
			else
			{
				uint32_t inv_light = 256 - light;
				uint32_t inv_desaturate = 256 - constants.desaturate;

				uint32_t intensity = ((red * 77 + green * 143 + blue * 37) >> 8) * constants.desaturate;

				red = (red * inv_desaturate + intensity) / 256;
				green = (green * inv_desaturate + intensity) / 256;
				blue = (blue * inv_desaturate + intensity) / 256;

				red = (constants.fade_red * inv_light + red * light) / 256;
				green = (constants.fade_green * inv_light + green * light) / 256;
				blue = (constants.fade_blue * inv_light + blue * light) / 256;

				red = (red * constants.light_red) / 256;
				green = (green * constants.light_green) / 256;
				blue = (blue * constants.light_blue) / 256;
			}

			dest[x] = 0xff000000 | (red << 16) | (green << 8) | blue;
		} while (++x <= x2);
	}
};

ApplySpecialColormapRGBACommand::ApplySpecialColormapRGBACommand(FSpecialColormap *colormap, DFrameBuffer *screen)
{
	buffer = screen->GetBuffer();
	pitch = screen->GetPitch();
	width = screen->GetWidth();
	height = screen->GetHeight();

	start_red = (int)(colormap->ColorizeStart[0] * 255);
	start_green = (int)(colormap->ColorizeStart[1] * 255);
	start_blue = (int)(colormap->ColorizeStart[2] * 255);
	end_red = (int)(colormap->ColorizeEnd[0] * 255);
	end_green = (int)(colormap->ColorizeEnd[1] * 255);
	end_blue = (int)(colormap->ColorizeEnd[2] * 255);
}

#ifdef NO_SSE
void ApplySpecialColormapRGBACommand::Execute(DrawerThread *thread)
{
	int y = thread->skipped_by_thread(0);
	int count = thread->count_for_thread(0, height);
	while (count > 0)
	{
		BYTE *pixels = buffer + y * pitch * 4;
		for (int x = 0; x < width; x++)
		{
			int fg_red = pixels[2];
			int fg_green = pixels[1];
			int fg_blue = pixels[0];

			int gray = (fg_red * 77 + fg_green * 143 + fg_blue * 37) >> 8;
			gray += (gray >> 7); // gray*=256/255
			int inv_gray = 256 - gray;

			int red = clamp((start_red * inv_gray + end_red * gray) >> 8, 0, 255);
			int green = clamp((start_green * inv_gray + end_green * gray) >> 8, 0, 255);
			int blue = clamp((start_blue * inv_gray + end_blue * gray) >> 8, 0, 255);

			pixels[0] = (BYTE)blue;
			pixels[1] = (BYTE)green;
			pixels[2] = (BYTE)red;
			pixels[3] = 0xff;

			pixels += 4;
		}
		y += thread->num_cores;
		count--;
	}
}
#else
void ApplySpecialColormapRGBACommand::Execute(DrawerThread *thread)
{
	int y = thread->skipped_by_thread(0);
	int count = thread->count_for_thread(0, height);
	__m128i gray_weight = _mm_set_epi16(256, 77, 143, 37, 256, 77, 143, 37);
	__m128i start_end = _mm_set_epi16(255, start_red, start_green, start_blue, 255, end_red, end_green, end_blue);
	while (count > 0)
	{
		BYTE *pixels = buffer + y * pitch * 4;
		int sse_length = width / 4;
		for (int x = 0; x < sse_length; x++)
		{
			// Unpack to integers:
			__m128i p = _mm_loadu_si128((const __m128i*)pixels);

			__m128i p16_0 = _mm_unpacklo_epi8(p, _mm_setzero_si128());
			__m128i p16_1 = _mm_unpackhi_epi8(p, _mm_setzero_si128());

			// Add gray weighting to colors
			__m128i mullo0 = _mm_mullo_epi16(p16_0, gray_weight);
			__m128i mullo1 = _mm_mullo_epi16(p16_1, gray_weight);
			__m128i p32_0 = _mm_unpacklo_epi16(mullo0, _mm_setzero_si128());
			__m128i p32_1 = _mm_unpackhi_epi16(mullo0, _mm_setzero_si128());
			__m128i p32_2 = _mm_unpacklo_epi16(mullo1, _mm_setzero_si128());
			__m128i p32_3 = _mm_unpackhi_epi16(mullo1, _mm_setzero_si128());

			// Transpose to get color components in individual vectors:
			__m128 tmpx = _mm_castsi128_ps(p32_0);
			__m128 tmpy = _mm_castsi128_ps(p32_1);
			__m128 tmpz = _mm_castsi128_ps(p32_2);
			__m128 tmpw = _mm_castsi128_ps(p32_3);
			_MM_TRANSPOSE4_PS(tmpx, tmpy, tmpz, tmpw);
			__m128i blue = _mm_castps_si128(tmpx);
			__m128i green = _mm_castps_si128(tmpy);
			__m128i red = _mm_castps_si128(tmpz);
			__m128i alpha = _mm_castps_si128(tmpw);

			// Calculate gray and 256-gray values:
			__m128i gray = _mm_srli_epi32(_mm_add_epi32(_mm_add_epi32(red, green), blue), 8);
			__m128i inv_gray = _mm_sub_epi32(_mm_set1_epi32(256), gray);

			// p32 = start * inv_gray + end * gray:
			__m128i gray0 = _mm_shuffle_epi32(gray, _MM_SHUFFLE(0, 0, 0, 0));
			__m128i gray1 = _mm_shuffle_epi32(gray, _MM_SHUFFLE(1, 1, 1, 1));
			__m128i gray2 = _mm_shuffle_epi32(gray, _MM_SHUFFLE(2, 2, 2, 2));
			__m128i gray3 = _mm_shuffle_epi32(gray, _MM_SHUFFLE(3, 3, 3, 3));
			__m128i inv_gray0 = _mm_shuffle_epi32(inv_gray, _MM_SHUFFLE(0, 0, 0, 0));
			__m128i inv_gray1 = _mm_shuffle_epi32(inv_gray, _MM_SHUFFLE(1, 1, 1, 1));
			__m128i inv_gray2 = _mm_shuffle_epi32(inv_gray, _MM_SHUFFLE(2, 2, 2, 2));
			__m128i inv_gray3 = _mm_shuffle_epi32(inv_gray, _MM_SHUFFLE(3, 3, 3, 3));
			__m128i gray16_0 = _mm_packs_epi32(gray0, inv_gray0);
			__m128i gray16_1 = _mm_packs_epi32(gray1, inv_gray1);
			__m128i gray16_2 = _mm_packs_epi32(gray2, inv_gray2);
			__m128i gray16_3 = _mm_packs_epi32(gray3, inv_gray3);
			__m128i gray16_0_mullo = _mm_mullo_epi16(gray16_0, start_end);
			__m128i gray16_1_mullo = _mm_mullo_epi16(gray16_1, start_end);
			__m128i gray16_2_mullo = _mm_mullo_epi16(gray16_2, start_end);
			__m128i gray16_3_mullo = _mm_mullo_epi16(gray16_3, start_end);
			__m128i gray16_0_mulhi = _mm_mulhi_epi16(gray16_0, start_end);
			__m128i gray16_1_mulhi = _mm_mulhi_epi16(gray16_1, start_end);
			__m128i gray16_2_mulhi = _mm_mulhi_epi16(gray16_2, start_end);
			__m128i gray16_3_mulhi = _mm_mulhi_epi16(gray16_3, start_end);
			p32_0 = _mm_srli_epi32(_mm_add_epi32(_mm_unpacklo_epi16(gray16_0_mullo, gray16_0_mulhi), _mm_unpackhi_epi16(gray16_0_mullo, gray16_0_mulhi)), 8);
			p32_1 = _mm_srli_epi32(_mm_add_epi32(_mm_unpacklo_epi16(gray16_1_mullo, gray16_1_mulhi), _mm_unpackhi_epi16(gray16_1_mullo, gray16_1_mulhi)), 8);
			p32_2 = _mm_srli_epi32(_mm_add_epi32(_mm_unpacklo_epi16(gray16_2_mullo, gray16_2_mulhi), _mm_unpackhi_epi16(gray16_2_mullo, gray16_2_mulhi)), 8);
			p32_3 = _mm_srli_epi32(_mm_add_epi32(_mm_unpacklo_epi16(gray16_3_mullo, gray16_3_mulhi), _mm_unpackhi_epi16(gray16_3_mullo, gray16_3_mulhi)), 8);

			p16_0 = _mm_packs_epi32(p32_0, p32_1);
			p16_1 = _mm_packs_epi32(p32_2, p32_3);
			p = _mm_packus_epi16(p16_0, p16_1);

			_mm_storeu_si128((__m128i*)pixels, p);
			pixels += 16;
		}

		for (int x = sse_length * 4; x < width; x++)
		{
			int fg_red = pixels[2];
			int fg_green = pixels[1];
			int fg_blue = pixels[0];

			int gray = (fg_red * 77 + fg_green * 143 + fg_blue * 37) >> 8;
			gray += (gray >> 7); // gray*=256/255
			int inv_gray = 256 - gray;

			int red = clamp((start_red * inv_gray + end_red * gray) >> 8, 0, 255);
			int green = clamp((start_green * inv_gray + end_green * gray) >> 8, 0, 255);
			int blue = clamp((start_blue * inv_gray + end_blue * gray) >> 8, 0, 255);

			pixels[0] = (BYTE)blue;
			pixels[1] = (BYTE)green;
			pixels[2] = (BYTE)red;
			pixels[3] = 0xff;

			pixels += 4;
		}

		y += thread->num_cores;
		count--;
	}
}
#endif

/////////////////////////////////////////////////////////////////////////////

void R_BeginDrawerCommands()
{
	DrawerCommandQueue::Begin();
}

void R_EndDrawerCommands()
{
	DrawerCommandQueue::End();
}

void R_DrawColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawColumnRGBACommand>();
}

void R_FillColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<FillColumnRGBACommand>();
}

void R_FillAddColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<FillAddColumnRGBACommand>();
}

void R_FillAddClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<FillAddClampColumnRGBACommand>();
}

void R_FillSubClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<FillSubClampColumnRGBACommand>();
}

void R_FillRevSubClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<FillRevSubClampColumnRGBACommand>();
}

void R_DrawFuzzColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawFuzzColumnRGBACommand>();
	fuzzpos = (fuzzpos + dc_yh - dc_yl) % FUZZTABLE;
}

void R_DrawAddColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawAddColumnRGBACommand>();
}

void R_DrawTranslatedColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawTranslatedColumnRGBACommand>();
}

void R_DrawTlatedAddColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawTlatedAddColumnRGBACommand>();
}

void R_DrawShadedColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawShadedColumnRGBACommand>();
}

void R_DrawAddClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawAddClampColumnRGBACommand>();
}

void R_DrawAddClampTranslatedColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawAddClampTranslatedColumnRGBACommand>();
}

void R_DrawSubClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSubClampColumnRGBACommand>();
}

void R_DrawSubClampTranslatedColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSubClampTranslatedColumnRGBACommand>();
}

void R_DrawRevSubClampColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawRevSubClampColumnRGBACommand>();
}

void R_DrawRevSubClampTranslatedColumn_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawRevSubClampTranslatedColumnRGBACommand>();
}

void R_DrawSpan_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanRGBACommand>();
}

void R_DrawSpanMasked_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanMaskedRGBACommand>();
}

void R_DrawSpanTranslucent_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanTranslucentRGBACommand>();
}

void R_DrawSpanMaskedTranslucent_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanMaskedTranslucentRGBACommand>();
}

void R_DrawSpanAddClamp_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanAddClampRGBACommand>();
}

void R_DrawSpanMaskedAddClamp_rgba()
{
	DrawerCommandQueue::QueueCommand<DrawSpanMaskedAddClampRGBACommand>();
}

void R_FillSpan_rgba()
{
	DrawerCommandQueue::QueueCommand<FillSpanRGBACommand>();
}

//extern FTexture *rw_pic; // For the asserts below

DWORD vlinec1_rgba()
{
	/*DWORD fracstep = dc_iscale;
	DWORD frac = dc_texturefrac;
	DWORD height = rw_pic->GetHeight();
	assert((frac >> vlinebits) < height);
	frac += (dc_count-1) * fracstep;
	assert((frac >> vlinebits) <= height);*/

	DrawerCommandQueue::QueueCommand<Vlinec1RGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void vlinec4_rgba()
{
	DrawerCommandQueue::QueueCommand<Vlinec4RGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

DWORD mvlinec1_rgba()
{
	DrawerCommandQueue::QueueCommand<Mvlinec1RGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void mvlinec4_rgba()
{
	DrawerCommandQueue::QueueCommand<Mvlinec4RGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

fixed_t tmvline1_add_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline1AddRGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void tmvline4_add_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline4AddRGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

fixed_t tmvline1_addclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline1AddClampRGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void tmvline4_addclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline4AddClampRGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

fixed_t tmvline1_subclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline1SubClampRGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void tmvline4_subclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline4SubClampRGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

fixed_t tmvline1_revsubclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline1RevSubClampRGBACommand>();
	return dc_texturefrac + dc_count * dc_iscale;
}

void tmvline4_revsubclamp_rgba()
{
	DrawerCommandQueue::QueueCommand<Tmvline4RevSubClampRGBACommand>();
	for (int i = 0; i < 4; i++)
		vplce[i] += vince[i] * dc_count;
}

void R_DrawFogBoundarySection_rgba(int y, int y2, int x1)
{
	for (; y < y2; ++y)
	{
		int x2 = spanend[y];
		DrawerCommandQueue::QueueCommand<DrawFogBoundaryLineRGBACommand>(y, x1, x2);
	}
}

void R_DrawFogBoundary_rgba(int x1, int x2, short *uclip, short *dclip)
{
	// To do: we do not need to create new spans when using rgba output - instead we should calculate light on a per pixel basis

	// This is essentially the same as R_MapVisPlane but with an extra step
	// to create new horizontal spans whenever the light changes enough that
	// we need to use a new colormap.

	double lightstep = rw_lightstep;
	double light = rw_light + rw_lightstep*(x2 - x1 - 1);
	int x = x2 - 1;
	int t2 = uclip[x];
	int b2 = dclip[x];
	int rcolormap = GETPALOOKUP(light, wallshade);
	int lcolormap;
	BYTE *basecolormapdata = basecolormap->Maps;

	if (b2 > t2)
	{
		clearbufshort(spanend + t2, b2 - t2, x);
	}

	R_SetColorMapLight(basecolormap, (float)light, wallshade);

	BYTE *fake_dc_colormap = basecolormap->Maps + (GETPALOOKUP(light, wallshade) << COLORMAPSHIFT);

	for (--x; x >= x1; --x)
	{
		int t1 = uclip[x];
		int b1 = dclip[x];
		const int xr = x + 1;
		int stop;

		light -= rw_lightstep;
		lcolormap = GETPALOOKUP(light, wallshade);
		if (lcolormap != rcolormap)
		{
			if (t2 < b2 && rcolormap != 0)
			{ // Colormap 0 is always the identity map, so rendering it is
			  // just a waste of time.
				R_DrawFogBoundarySection_rgba(t2, b2, xr);
			}
			if (t1 < t2) t2 = t1;
			if (b1 > b2) b2 = b1;
			if (t2 < b2)
			{
				clearbufshort(spanend + t2, b2 - t2, x);
			}
			rcolormap = lcolormap;
			R_SetColorMapLight(basecolormap, (float)light, wallshade);
			fake_dc_colormap = basecolormap->Maps + (GETPALOOKUP(light, wallshade) << COLORMAPSHIFT);
		}
		else
		{
			if (fake_dc_colormap != basecolormapdata)
			{
				stop = MIN(t1, b2);
				while (t2 < stop)
				{
					int y = t2++;
					DrawerCommandQueue::QueueCommand<DrawFogBoundaryLineRGBACommand>(y, xr, spanend[y]);
				}
				stop = MAX(b1, t2);
				while (b2 > stop)
				{
					int y = --b2;
					DrawerCommandQueue::QueueCommand<DrawFogBoundaryLineRGBACommand>(y, xr, spanend[y]);
				}
			}
			else
			{
				t2 = MAX(t2, MIN(t1, b2));
				b2 = MIN(b2, MAX(b1, t2));
			}

			stop = MIN(t2, b1);
			while (t1 < stop)
			{
				spanend[t1++] = x;
			}
			stop = MAX(b2, t2);
			while (b1 > stop)
			{
				spanend[--b1] = x;
			}
		}

		t2 = uclip[x];
		b2 = dclip[x];
	}
	if (t2 < b2 && rcolormap != 0)
	{
		R_DrawFogBoundarySection_rgba(t2, b2, x1);
	}
}
