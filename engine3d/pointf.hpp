#ifndef OCTILLION_ENGINE3D_POINTF_HEADER
#define OCTILLION_ENGINE3D_POINTF_HEADER

#include <cmath>

#include "matrix.hpp"

class PointF
{
public:
	PointF() {};

	PointF(double x, double y, double z) : x_(x), y_(y), z_(z) {}

  PointF(double x, double y, double z, int_fast32_t u, int_fast32_t v) : x_(x), y_(y), z_(z), u_(u), v_(v) {}

	PointF(const PointF& rhs) : x_(rhs.x_), y_(rhs.y_), z_(rhs.z_), u_(rhs.u_), v_(rhs.v_) {}

	PointF& operator=(const PointF& rhs)
	{
		x_ = rhs.x_;
		y_ = rhs.y_;
		z_ = rhs.z_;
    u_ = rhs.u_;
    v_ = rhs.v_;
		return *this;
	}

	void set(double x, double y, double z)
	{
		x_ = x;
		y_ = y;
		z_ = z;
	}

  void setuv(int_fast32_t u, int_fast32_t v)
  {
    u_ = u;
    v_ = v;
  }

	double x() const { return x_; }
	double y() const { return y_; }
	double z() const { return z_; }

	void x(double x) { x_ = x; }
	void y(double y) { y_ = y; }
	void z(double z) { z_ = z; }

	friend PointF operator - (const PointF& lhs, const PointF& rhs)
	{
		PointF ret(lhs.x_ - rhs.x_, lhs.y_ - rhs.y_, lhs.z_ - rhs.z_);
		return ret;
	}

	friend PointF operator * (const PointF& lhs, const Matrix<double>& matrix)
	{
		Matrix<double> lhs_matrix(lhs.x_, lhs.y_, lhs.z_);
		Matrix<double> result = lhs_matrix * matrix;
		PointF pt(result.at(0), result.at(1), result.at(2), lhs.u_, lhs.v_);
		return pt;
	}

public:
	double x_, y_, z_;

  // UV coordinate
  int_fast32_t u_, v_;
};


#endif