#pragma once

#include "resource.h"

class ScreenBuffer
{
public:
	//Singleton
	static ScreenBuffer& get_instance()
	{
		static ScreenBuffer instance;
		return instance;
	}

	// avoid accidentally copy
	ScreenBuffer(ScreenBuffer const&) = delete;
	void operator = (ScreenBuffer const&) = delete;

private:
	ScreenBuffer() {}
	~ScreenBuffer()
	{
		if (buffer_ != NULL)
		{
			delete[] buffer_;
		}
	}

public:
	BYTE* get(LONG width, LONG height, WORD bytesPerPixel)
	{
		if (width != width_ || height != height_)
		{
			width_ = width;
			height_ = height;

			if (buffer_ != NULL)
			{
				delete[] buffer_;
			}

			buffer_ = new BYTE[width_ * height_ * bytesPerPixel_];
		}

		return buffer_;
	}

private:
	LONG width_ = 0, height_ = 0;
	WORD  bytesPerPixel_ = 4; // default 4 bytes
	BYTE* buffer_ = NULL;
};