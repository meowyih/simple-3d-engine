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
	void paintCubePoints(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);
	void paintCubeLines(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

private:
	Sample_1();
	~Sample_1() {}

private:
	unsigned int rotate_degree_ = 0;
	PointF pts_[8];
};

class Sample_2
{
private:
	double pi_ = 3.1415926;

public:

	//Singleton
	static Sample_2& get_instance()
	{
		static Sample_2 instance;
		return instance;
	}

	// avoid accidentally copy
	Sample_2(Sample_2 const&) = delete;
	void operator = (Sample_2 const&) = delete;

	// paint
	void paintFace(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

private:
	Sample_2();
	~Sample_2() {}

private:
	unsigned int rotate_degree_ = 0;
	Face face_;
};