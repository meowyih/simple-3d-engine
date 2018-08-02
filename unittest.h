#pragma once

#include "engine3d/pointf.hpp"
#include "engine3d/face.hpp"

class Sample_1
{
private:
	double pi_ = 3.1415926;

public:

	//Singleton
	static Sample_1& get_instance()
	{
		static Sample_1 instance;
		return instance;
	}

	// avoid accidentally copy
	Sample_1(Sample_1 const&) = delete;
	void operator = (Sample_1 const&) = delete;

	// paint
	void paint(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

private:
	Sample_1();
	~Sample_1() {}

private:
	unsigned int rotate_degree_ = 0;
	PointF pts_[8];
};