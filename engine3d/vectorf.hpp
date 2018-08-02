#ifndef OCTILLION_VECTOR_HEADER
#define OCTILLION_VECTOR_HEADER

#include <cmath>

#include "pointf.hpp"
#include "matrix.hpp"

class VectorF : public PointF
{
public:
	// default constructor
	VectorF() {};

	// constructor by x, y and z
	VectorF(double x, double y, double z) : PointF(x, y, z) {}

	// construct by PointF
	VectorF(const PointF& pf) : PointF(pf) {}

	// copy constructor
	VectorF(const VectorF& rhs)
	{
		x_ = rhs.x_;
		y_ = rhs.y_;
		z_ = rhs.z_;
	}

	// copy operator
	VectorF& operator=(const VectorF& rhs)
	{
		x_ = rhs.x_;
		y_ = rhs.y_;
		z_ = rhs.z_;
		return *this;
	}

	// translate the vector by matrix
	friend VectorF operator* (const VectorF& lhs, const Matrix<double>& rhs)
	{
		Matrix<double> lhs_matrix(lhs.x_, lhs.y_, lhs.z_);
		Matrix<double> result = lhs_matrix * rhs;
		return VectorF(result.at(0), result.at(1), result.at(2));
	}

	// get the normalized vector
	static VectorF normalize(const VectorF& vector)
	{
		double magnitude = sqrt(vector.x_ * vector.x_ + vector.y_ * vector.y_ + vector.z_*vector.z_);
		return VectorF(vector.x_ / magnitude, vector.y_ / magnitude, vector.z_ / magnitude);

	}

	// calculate the length of a vector
	static double magnitude(const VectorF& vector)
	{
		return sqrt(vector.x_ * vector.x_ + vector.y_ * vector.y_ + vector.z_*vector.z_);
	}

	// calculate cosine between v1 and v2_ by applying dot product
	static double cosine(VectorF v1, VectorF v2)
	{
		double cosine = (v1.x_ * v2.x_ + v1.y_ * v2.y_ + v1.z_ * v2.z_);

		return cosine;
	}
};

#endif