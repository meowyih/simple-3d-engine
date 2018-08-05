#ifndef OCTILLION_ENGINE3D_POINT_HEADER
#define OCTILLION_ENGINE3D_POINT_HEADER

#include <cstdint>

#include "pointf.hpp"
#include "vectorf.hpp"
#include "matrix.hpp"

class PointI
{
public:
  static const int_fast32_t COS_PRECISENESS = 100000;
public:
	PointI() {};
	PointI(int_fast32_t x, int_fast32_t y, int_fast32_t z) : x_(x), y_(y), z_(z) {}
  PointI(int_fast32_t x, int_fast32_t y, int_fast32_t z, int_fast32_t cos ) : x_(x), y_(y), z_(z), cos_(cos) {}

  PointI(const PointI& rhs) : x_(rhs.x_), y_(rhs.y_), z_(rhs.z_), cos_(rhs.cos_), u_(rhs.u_), v_(rhs.v_) {}

	PointI(const PointF& ptf)
	{
		x_ = (int_fast32_t)round(ptf.x());
		y_ = (int_fast32_t)round(ptf.y());
		z_ = (int_fast32_t)round(ptf.z());
    u_ = ptf.u_;
    v_ = ptf.v_;
	}

	int_fast32_t x() const { return x_; }
	int_fast32_t y() const { return y_; }
	int_fast32_t z() const { return z_; }

	void set(int_fast32_t x, int_fast32_t y, int_fast32_t z)
	{
		x_ = x;
		y_ = y;
		z_ = z;
	}

	void x(int_fast32_t x) { x_ = x; }
	void y(int_fast32_t y) { y_ = y; }
	void z(int_fast32_t z) { z_ = z; }

	friend PointI operator - (const PointI& lhs, const PointI& rhs)
	{
		PointI ret(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_);
		return ret;
	}

public:
	int_fast32_t x_, y_, z_;

  // Phong shading usage, cosine value times 1000
  // this value should between -1000 to 1000
  int_fast32_t cos_; 

  // UV coordinate
  int_fast32_t u_, v_;
};


#endif