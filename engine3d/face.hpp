#ifndef OCTILLION_ENGINE3D_FACE_HEADER
#define OCTILLION_ENGINE3D_FACE_HEADER

#include <vector>

#include "point.hpp"
#include "pointf.hpp"
#include "vectorf.hpp"
#include "utility.hpp"
#include "matrix.hpp"

class Face
{
public:
	Face() {}

	Face(PointF pt1, PointF pt2, PointF pt3) : pt1_(pt1), pt2_(pt2), pt3_(pt3) {}

	Face(PointF pt1, PointF pt1n,
		PointF pt2, PointF pt2n,
		PointF pt3, PointF pt3n) :
		pt1_(pt1), pt2_(pt2), pt3_(pt3),
		pt1n_(VectorF::normalize(pt1n)),
		pt2n_(VectorF::normalize(pt2n)),
		pt3n_(VectorF::normalize(pt3n)) {}

	void set(PointF pt1, PointF pt2, PointF pt3)
	{
		points_.clear();
		pt1_ = pt1;
		pt2_ = pt2;
		pt3_ = pt3;
	}

	void rasterize()
	{
		points_.clear();
		Utility::rasterize_triangle(
			static_cast<PointI>(pt1_), 
			static_cast<PointI>(pt2_), 
			static_cast<PointI>(pt3_), 
			points_);
	}

	void rasterize(PointF light, std::vector<double>& cosine)
	{
		points_.clear();
		cosine.clear();
		Utility::rasterize_triangle(
			light,
			pt1_, pt1n_,
			pt2_, pt2n_,
			pt3_, pt3n_,
			points_,
			cosine );
		return;
	}

	friend Face operator* (const Face& lhs, const Matrix<double>& rhs)
	{
		Matrix<double> pt1_matrix(lhs.pt1_.x_, lhs.pt1_.y_, lhs.pt1_.z_);
		Matrix<double> pt1_trans_matrix = pt1_matrix * rhs;
		Matrix<double> pt2_matrix(lhs.pt2_.x_, lhs.pt2_.y_, lhs.pt2_.z_);
		Matrix<double> pt2_trans_matrix = pt2_matrix * rhs;
		Matrix<double> pt3_matrix(lhs.pt3_.x_, lhs.pt3_.y_, lhs.pt3_.z_);
		Matrix<double> pt3_trans_matrix = pt3_matrix * rhs;

		Matrix<double> pt1n_matrix(lhs.pt1n_.x_, lhs.pt1n_.y_, lhs.pt1n_.z_);
		Matrix<double> pt1n_trans_matrix = pt1n_matrix * rhs;
		Matrix<double> pt2n_matrix(lhs.pt2n_.x_, lhs.pt2n_.y_, lhs.pt2n_.z_);
		Matrix<double> pt2n_trans_matrix = pt2n_matrix * rhs;
		Matrix<double> pt3n_matrix(lhs.pt3n_.x_, lhs.pt3n_.y_, lhs.pt3n_.z_);
		Matrix<double> pt3n_trans_matrix = pt3n_matrix * rhs;

		return Face(
			PointF(pt1_trans_matrix.at(0), pt1_trans_matrix.at(1), pt1_trans_matrix.at(2)),
			VectorF(pt1n_trans_matrix.at(0), pt1n_trans_matrix.at(1), pt1n_trans_matrix.at(2)),
			PointF(pt2_trans_matrix.at(0), pt2_trans_matrix.at(1), pt2_trans_matrix.at(2)),
			VectorF(pt2n_trans_matrix.at(0), pt2n_trans_matrix.at(1), pt2n_trans_matrix.at(2)),
			PointF(pt3_trans_matrix.at(0), pt3_trans_matrix.at(1), pt3_trans_matrix.at(2)),
			VectorF(pt3n_trans_matrix.at(0), pt3n_trans_matrix.at(1), pt3n_trans_matrix.at(2))
		);

	}

public:
	PointF pt1_, pt2_, pt3_;
	VectorF pt1n_, pt2n_, pt3n_;

	std::vector<PointI> points_;
};

#endif