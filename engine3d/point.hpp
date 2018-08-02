#ifndef OCTILLION_ENGINE3D_POINT_HEADER
#define OCTILLION_ENGINE3D_POINT_HEADER

#include <cstdint>

#include "pointf.hpp"
#include "matrix.hpp"


class Point
{
public:
	Point() {};
	Point(int_fast32_t x, int_fast32_t y, int_fast32_t z) : x_(x), y_(y), z_(z) {}

	Point(const PointF& ptf)
	{
		x_ = (int_fast32_t)round(ptf.x());
		y_ = (int_fast32_t)round(ptf.y());
		z_ = (int_fast32_t)round(ptf.z());
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

	friend Point operator - (const Point& lhs, const Point& rhs)
	{
		Point ret(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_);
		return ret;
	}

public:
	int_fast32_t x_, y_, z_;
};


#endif