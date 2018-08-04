#pragma once

#include "engine3d/pointf.hpp"
#include "engine3d/face.hpp"
#include "engine3d/mesh.hpp"

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

class Sample_3
{
private:
  double pi_ = 3.1415926;

public:

  //Singleton
  static Sample_3& get_instance()
  {
    static Sample_3 instance;
    return instance;
  }

  // avoid accidentally copy
  Sample_3(Sample_3 const&) = delete;
  void operator = (Sample_3 const&) = delete;

  // paint
  void paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

private:
  Sample_3();
  ~Sample_3() {}

private:
  unsigned int rotate_degree_ = 0;
  Mesh mesh_;
};

class Sample_4
{
private:
  double pi_ = 3.1415926;

public:

  //Singleton
  static Sample_4& get_instance()
  {
    static Sample_4 instance;
    return instance;
  }

  // avoid accidentally copy
  Sample_4(Sample_4 const&) = delete;
  void operator = (Sample_4 const&) = delete;

  // paint
  void paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

  // help static function to get the center point for a face
  static PointF center(Face face);

private:
  Sample_4();
  ~Sample_4() {}

private:
  unsigned int rotate_degree_ = 0;
  Face abc_, abd_, acd_, bcd_;
  VectorF n_abc_, n_abd_, n_acd_, n_bcd_;
};

class Sample_5
{
private:
  double pi_ = 3.1415926;

public:

  //Singleton
  static Sample_5& get_instance()
  {
    static Sample_5 instance;
    return instance;
  }

  // avoid accidentally copy
  Sample_5(Sample_5 const&) = delete;
  void operator = (Sample_5 const&) = delete;

  // paint
  void paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel);

private:
  Sample_5();
  ~Sample_5() {}

private:
  unsigned int rotate_degree_ = 0;
  Mesh mesh_;
};