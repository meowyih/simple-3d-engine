#ifndef OCTILLION_UTILITY_HEADER
#define OCTILLION_UTILITY_HEADER

#include <vector>
#include <cstdint>   // int_fast32_t
#include <cmath>     // abs()
#include <algorithm> // max()

#include "matrix.hpp"
#include "point.hpp"
#include "pointf.hpp"
#include "vectorf.hpp"

class Utility
{
public:

	static void reset_matrix(Matrix<double>& matrix)
	{
		double data[] = { 1.0, 0, 0, 0,
			0, 1.0, 0, 0,
			0, 0, 1.0, 0,
			0, 0, 0, 1.0 };

		matrix.set(4, 4, data);
	}

	static void set_scale(Matrix<double>& matrix, double scale)
	{
		double data[16] = {
			scale, 0, 0, 0,
			0, scale, 0, 0,
			0, 0, scale, 0,
			0, 0, 0,     1 };

		matrix.set(data);
	}

	static void set_translate(Matrix<double>& matrix, double x, double y, double z)
	{
		double data[16] = {
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			x, y, z, 1 };

		matrix.set(data);
	}

	static void set_rotateX(Matrix<double>& matrix, double radian)
	{
		double data[16] = {
			1, 0, 0, 0,
			0, cos(radian), sin(radian), 0,
			0, -sin(radian), cos(radian), 0,
			0, 0, 0, 1 };

		matrix.set(data);
	}

	static void set_rotateY(Matrix<double>& matrix, double radian)
	{
		double data[16] = {
			cos(radian), 0, -sin(radian), 0,
			0, 1, 0, 0,
			sin(radian), 0, cos(radian), 0,
			0, 0, 0, 1 };

		matrix.set(data);
	}

	static void set_rotateZ(Matrix<double>& matrix, double radian)
	{
		double data[16] = {
			cos(radian), sin(radian), 0, 0,
			-sin(radian), cos(radian), 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1 };

		matrix.set(data);
	}

	// standard bresenham 2d
	inline static void bresenham2d(Point pt1, Point pt2, std::vector<Point>& out)
	{
		int_fast32_t x0 = pt1.x_;
		int_fast32_t y0 = pt1.y_;
		int_fast32_t x1 = pt2.x_;
		int_fast32_t y1 = pt2.y_;
		int_fast32_t tmp;
		bool steep = abs(pt2.y_ - pt1.y_) > abs(pt2.x_ - pt1.x_) ? true : false;
		Point from, to;
		if (steep)
		{
			tmp = x0;
			x0 = y0;
			y0 = tmp;

			tmp = x1;
			x1 = y1;
			y1 = tmp;
		}

		if (x0 > x1)
		{
			tmp = x0;
			x0 = x1;
			x1 = tmp;

			tmp = y0;
			y0 = y1;
			y1 = tmp;
		}

		int_fast32_t deltax = x1 - x0;
		int_fast32_t deltay = abs(y1 - y0);
		int_fast32_t error = deltax >> 1;
		int_fast32_t ystep;
		int_fast32_t y = y0;

		if (y0 < y1)
			ystep = 1;
		else
			ystep = -1;

		for (int_fast32_t x = x0; x <= x1; x++)
		{
			if (steep)
			{
				Point pt(y, x, pt1.z_);
				out.push_back(pt);
			}
			else
			{
				Point pt(x, y, pt1.z_);
				out.push_back(pt);
			}

			error -= deltay;

			if (error < 0)
			{
				y += ystep;
				error += deltax;
			}
		}
	}

	inline static void bresenham2d2(Point pt1, Point pt2, std::vector<Point>& out)
	{
		Point pos(pt1);
		Point delta = pt2 - pt1;
		Point abs_delta(abs(delta.x_), abs(delta.y_), abs(delta.z_));
		Point double_delta(abs_delta.x_ << 1, abs_delta.y_ << 1, abs_delta.z_ << 1);

		int_fast32_t x_inc, y_inc, z_inc;

		out.reserve(abs_delta.x_ + abs_delta.y_);

		x_inc = (delta.x_ < 0) ? -1 : 1;
		y_inc = (delta.y_ < 0) ? -1 : 1;
		z_inc = (delta.z_ < 0) ? -1 : 1;

		if (abs_delta.x_ >= abs_delta.y_)
		{
			int_fast32_t err_1 = double_delta.y_ - abs_delta.x_;
			int_fast32_t err_2 = double_delta.z_ - abs_delta.x_;

			for (int_fast32_t i = 0; i < abs_delta.x_; i++)
			{
				out.push_back(pos);

				if (err_1 > 0)
				{
					pos.y_ += y_inc;
					err_1 -= double_delta.x_;
				}

				if (err_2 > 0)
				{
					pos.z_ += z_inc;
					err_2 -= double_delta.x_;
				}

				err_1 += double_delta.y_;
				err_2 += double_delta.z_;

				pos.x_ += x_inc;
			}
		}
		else // if (abs_delta.y_ >= abs_delta.x_)
		{
			int_fast32_t err_1 = double_delta.x_ - abs_delta.y_;
			int_fast32_t err_2 = double_delta.z_ - abs_delta.y_;

			for (int_fast32_t i = 0; i < abs_delta.y_; i++)
			{
				out.push_back(pos);

				if (err_1 > 0)
				{
					pos.x_ += x_inc;
					err_1 -= double_delta.y_;
				}

				if (err_2 > 0)
				{
					pos.z_ += z_inc;
					err_2 -= double_delta.y_;
				}

				err_1 += double_delta.x_;
				err_2 += double_delta.z_;

				pos.y_ += y_inc;
			}
		}

		out.push_back(pos);
	}

	// line3d uses Bresenham's algorithm to generate the 3 dimensional points on a
	// line from (x1, y1, z1) to (x2, y2, z2), reference implementation found here
	// http://www.ict.griffith.edu.au/anthony/info/graphics/bresenham.procs (3D)
	// https://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm (2D)
	inline static void bresenham3d(Point pt1, Point pt2, std::vector<Point>& out)
	{
		Point pos(pt1);
		Point delta = pt2 - pt1;
		Point abs_delta(abs(delta.x_), abs(delta.y_), abs(delta.z_));
		Point double_delta(abs_delta.x_ << 1, abs_delta.y_ << 1, abs_delta.z_ << 1);

		int_fast32_t x_inc, y_inc, z_inc;

		out.reserve(abs_delta.x_ + abs_delta.y_ + abs_delta.z_);

		x_inc = (delta.x_ < 0) ? -1 : 1;
		y_inc = (delta.y_ < 0) ? -1 : 1;
		z_inc = (delta.z_ < 0) ? -1 : 1;

		if (abs_delta.x_ >= abs_delta.y_ && abs_delta.x_ >= abs_delta.z_)
		{
			int_fast32_t err_1 = double_delta.y_ - abs_delta.x_;
			int_fast32_t err_2 = double_delta.z_ - abs_delta.x_;

			for (int_fast32_t i = 0; i < abs_delta.x_; i++)
			{
				out.push_back(pos);

				if (err_1 > 0)
				{
					pos.y_ += y_inc;
					err_1 -= double_delta.x_;
				}

				if (err_2 > 0)
				{
					pos.z_ += z_inc;
					err_2 -= double_delta.x_;
				}

				err_1 += double_delta.y_;
				err_2 += double_delta.z_;

				pos.x_ += x_inc;
			}
		}
		else if (abs_delta.y_ >= abs_delta.x_ && abs_delta.y_ >= abs_delta.z_)
		{
			int_fast32_t err_1 = double_delta.x_ - abs_delta.y_;
			int_fast32_t err_2 = double_delta.z_ - abs_delta.y_;

			for (int_fast32_t i = 0; i < abs_delta.y_; i++)
			{
				out.push_back(pos);

				if (err_1 > 0)
				{
					pos.x_ += x_inc;
					err_1 -= double_delta.y_;
				}

				if (err_2 > 0)
				{
					pos.z_ += z_inc;
					err_2 -= double_delta.y_;
				}

				err_1 += double_delta.x_;
				err_2 += double_delta.z_;

				pos.y_ += y_inc;
			}
		}
		else
		{
			int_fast32_t err_1 = double_delta.y_ - abs_delta.z_;
			int_fast32_t err_2 = double_delta.x_ - abs_delta.z_;

			for (int_fast32_t i = 0; i < abs_delta.z_; i++)
			{
				// if x, y are the same, we don't have to store the pixel
				out.push_back(pos);
		
				if (err_1 > 0)
				{
					pos.y_ += y_inc;
					err_1 -= double_delta.z_;
				}

				if (err_2 > 0)
				{
					pos.x_ += x_inc;
					err_2 -= double_delta.z_;
				}

				err_1 += double_delta.y_;
				err_2 += double_delta.x_;

				pos.z_ += z_inc;
			}
		}

		out.push_back(pos);
	}

	inline static void rasterize_triangle(Point pt1, Point pt2, Point pt3, std::vector<Point>& out)
	{
		Point top, middle, bottom;

		std::vector<Point> vertices_top_middle;
		std::vector<Point> vertices_top_bottom;
		std::vector<Point> vertices_middle_bottom;

		if (pt1.x_ > pt2.x_ && pt1.x_ > pt3.x_ )
		{
			top = pt1;
			if (pt2.x_ > pt3.x_)
			{
				middle = pt2;
				bottom = pt3;
			}
			else
			{
				middle = pt3;
				bottom = pt2;
			}
		}
		else if (pt2.x_ > pt1.x_ && pt2.x_ > pt3.x_)
		{
			top = pt2;
			if (pt1.x_ > pt3.x_)
			{
				middle = pt1;
				bottom = pt3;
			}
			else
			{
				middle = pt3;
				bottom = pt1;
			}
		}
		else
		{
			top = pt3;
			if (pt1.x_ > pt2.x_)
			{
				middle = pt1;
				bottom = pt2;
			}
			else
			{
				middle = pt2;
				bottom = pt1;
			}
		}

		bresenham2d2(top, middle, vertices_top_middle);
		bresenham2d2(top, bottom, vertices_top_bottom);
		bresenham2d2(middle, bottom, vertices_middle_bottom);

		auto it1 = vertices_top_bottom.begin();
		auto it2 = vertices_top_middle.begin();
		auto it3 = vertices_middle_bottom.begin();

		int_fast32_t x_cur = it1->x_;

		size_t reserve = abs(top.x_ - bottom.x_ + 1) * (max(abs(top.y_ - middle.y_), abs(top.y_ - bottom.y_)) + 1);

		if (reserve > 4096)
		{
			out.reserve(4096);
		}
		else
		{
			out.reserve(reserve);
		}

		// fill upper triangle
		while (it2 != vertices_top_middle.end() && it1 != vertices_top_bottom.end())
		{
			if (it2->x_ == x_cur)
			{
				std::vector<Point> points;
				bresenham2d2(*it1, *it2, points);
				out.insert(out.end(), points.begin(), points.end());
				x_cur--;
			}

			// find next it1->x_ == it2->x_ == x_cur
			while (it1 != vertices_top_bottom.end())
			{					
				if (it1->x_ == x_cur)
					break;
				it1++;
			}

			while (it2 != vertices_top_middle.end())
			{					
				if (it2->x_ == x_cur)
					break;
				it2++;
			}
		}

		// fill lower triangle
		while (it3 != vertices_middle_bottom.end() && it1 != vertices_top_bottom.end())
		{
			if (it3->x_ == x_cur)
			{
				std::vector<Point> points;
				bresenham2d2(*it1, *it3, points);
				out.insert(out.end(), points.begin(), points.end());
				x_cur--;
			}

			// find next it1->x_ == it2->x_ == x_cur
			while (it1 != vertices_top_bottom.end())
			{					
				if (it1->x_ == x_cur)
					break;
				it1++;
			}

			while (it3 != vertices_middle_bottom.end())
			{					
				if (it3->x_ == x_cur)
					break;
				it3++;
			}
		}
	}

	inline static void rasterize_triangle(
		PointF light,
		PointF pt1, VectorF pt1n, 
		PointF pt2, VectorF pt2n,
		PointF pt3, VectorF pt3n,
		std::vector<Point>& out,
		std::vector<double>& cosine )
	{
		PointF top, middle, bottom;
		VectorF top_n, middle_n, bottom_n;
		double top_c, middle_c, bottom_c;
			
		std::vector<Point> vertices_top_middle;
		std::vector<Point> vertices_top_bottom;
		std::vector<Point> vertices_middle_bottom;

		if (pt1.x_ >= pt2.x_ && pt1.x_ >= pt3.x_)
		{
			top = pt1;
			top_n = pt1n;
			if (pt2.x_ > pt3.x_)
			{
				middle = pt2;
				middle_n = pt2n;
				bottom = pt3;
				bottom_n = pt3n;
			}
			else
			{
				middle = pt3;
				middle_n = pt3n;
				bottom = pt2;
				bottom_n = pt2n;
			}
		}
		else if (pt2.x_ >= pt1.x_ && pt2.x_ >= pt3.x_)
		{
			top = pt2;
			top_n = pt2n;
			if (pt1.x_ > pt3.x_)
			{
				middle = pt1;
				middle_n = pt1n;
				bottom = pt3;
				bottom_n = pt3n;
			}
			else
			{
				middle = pt3;
				middle_n = pt3n;
				bottom = pt1;
				bottom_n = pt1n;
			}
		}
		else
		{
			top = pt3;
			top_n = pt3n;
			if (pt1.x_ > pt2.x_)
			{
				middle = pt1;
				middle_n = pt1n;
				bottom = pt2;
				bottom_n = pt2n;
			}
			else
			{
				middle = pt2;
				middle_n = pt2n;
				bottom = pt1;
				bottom_n = pt1n;
			}
		}

		VectorF top_light = VectorF::normalize(light - top);
		VectorF middle_light = VectorF::normalize(light - middle);
		VectorF bottom_light = VectorF::normalize(light - bottom);

		top_c = VectorF::cosine(top_light, top_n);
		middle_c = VectorF::cosine(middle_light, middle_n);
		bottom_c = VectorF::cosine(bottom_light, bottom_n);

		// special case: 3 points has same x_
		if (top.x_ == middle.x_ && middle.x_ == bottom.x_)
		{
			std::vector<Point> vertices_y;
			PointF pt_ymin, pt_ymax;
			double pt_ymin_c, pt_ymax_c;

			if (top.y_ >= middle.y_ && top.y_ >= bottom.y_)
			{
				pt_ymax = top;
				pt_ymax_c = top_c;
			}
			else if (middle.y_ >= bottom.y_)
			{
				pt_ymax = middle;
				pt_ymax_c = middle_c;
			}
			else
			{
				pt_ymax = bottom;
				pt_ymax_c = bottom_c;
			}

			if (top.y_ <= middle.y_ && top.y_ <= bottom.y_)
			{
				pt_ymin = top;
				pt_ymin_c = top_c;
			}
			else if (middle.y_ <= bottom.y_)
			{
				pt_ymin = middle;
				pt_ymin_c = middle_c;
			}
			else
			{
				pt_ymin = bottom;
				pt_ymin_c = bottom_c;
			}

			Point pt_start(pt_ymin);
			Point pt_end(pt_ymax);
			bresenham2d2(pt_start, pt_end, vertices_y);

			for (Point pt : vertices_y)
			{
				double gradient1 = gradient(pt_ymin, pt_ymax, PointF(pt.x_, pt.y_, pt.z_));
				double cos = interpolation(pt_ymin_c, pt_ymax_c, gradient1);
				cosine.push_back(cos);
			}

			return;
		}

		// general cases
		bresenham2d2(static_cast<Point>(top), static_cast<Point>(middle), vertices_top_middle);
		bresenham2d2(static_cast<Point>(top), static_cast<Point>(bottom), vertices_top_bottom);
		bresenham2d2(static_cast<Point>(middle), static_cast<Point>(bottom), vertices_middle_bottom);

		auto it1 = vertices_top_bottom.begin();
		auto it2 = vertices_top_middle.begin();
		auto it3 = vertices_middle_bottom.begin();

		int_fast32_t x_cur = it1->x_;

		Point topi(top), bottomi(bottom), middlei(middle);
		size_t reserve = 
			abs((topi.x_ - bottomi.x_ + 1)
			* (max(abs(topi.y_ - middlei.y_) + 1, abs(topi.y_ - bottomi.y_))) + 1);

		if (reserve > 1024)
		{
			out.reserve(1024);
		}
		else
		{
			out.reserve(reserve);
		}

		// fill upper triangle, skip it if there is no upper triangle (single line)
		if (top.x_ != middle.x_)
		{
			while (it2 != vertices_top_middle.end() && it1 != vertices_top_bottom.end())
			{
				if (it2->x_ == x_cur)
				{
					std::vector<Point> points;

					double gradient1 = gradient(top, bottom, PointF(it1->x_, it1->y_, it1->z_));
					double gradient2 = gradient(top, middle, PointF(it2->x_, it2->y_, it2->z_));

					double start_c = interpolation(top_c, bottom_c, gradient1);
					double end_c = interpolation(top_c, middle_c, gradient2);

					bresenham2d2(*it1, *it2, points);

					for (Point pt : points)
					{
						double gradient =
							(pt.y_ == it1->y_ ? 1 : (it2->y_*1.0 - pt.y_*1.0) / (it2->y_*1.0 - it1->y_*1.0));

						if (gradient > 1)
						{
							gradient = 1;
						}

						double cos = interpolation(start_c, end_c, gradient);

						cosine.push_back(cos);
					}

					out.insert(out.end(), points.begin(), points.end());
					x_cur--;
				}

				// find next it2->x_ == x_cur
				while (it2 != vertices_top_middle.end())
				{
					if (it2->x_ == x_cur)
					{
						break;
					}
					it2++;
				}

				// if all it2 point are consumed, we retore x_cur and break;
				if (it2 == vertices_top_middle.end())
				{
					break;
				}
				else
				{
					// find next it1->x_ == x_cur
					while (it1 != vertices_top_bottom.end())
					{
						if (it1->x_ == x_cur)
							break;
						it1++;
					}
				}
			}
		}

		// fill lower triangle, skip it if there is no upper triangle (single line)
		if (middle.x_ != bottom.x_)
		{
			while (it3 != vertices_middle_bottom.end() && it1 != vertices_top_bottom.end())
			{
				if (it3->x_ == x_cur)
				{
					std::vector<Point> points;

					double gradient1 = gradient(top, bottom, PointF(it1->x_, it1->y_, it1->z_));
					double gradient3 = gradient(middle, bottom, PointF(it3->x_, it3->y_, it3->z_));

					double start_c = interpolation(top_c, bottom_c, gradient1);
					double end_c = interpolation(middle_c, bottom_c, gradient3);

					bresenham2d2(*it1, *it3, points);

					for (Point pt : points)
					{
						double gradient =
							(pt.y_ == it1->y_ ? 1 : (it3->y_*1.0 - pt.y_*1.0) / (it3->y_*1.0 - it1->y_*1.0));

						double cos = interpolation(start_c, end_c, gradient);
						cosine.push_back(cos);
					}

					out.insert(out.end(), points.begin(), points.end());
					x_cur--;
				}
				// find next it3->x_ == x_cur
				while (it3 != vertices_middle_bottom.end())
				{
					if (it3->x_ == x_cur)
						break;
					it3++;
				}

				if (it3 == vertices_middle_bottom.end())
				{
					// all it3 vertices were consumed
					break;
				}
				else
				{
					// find next it1->x_ == it2->x_ == x_cur
					while (it1 != vertices_top_bottom.end())
					{
						if (it1->x_ == x_cur)
							break;
						it1++;
					}
				}
			}
		}
	}
		
	static double interpolation(double v1, double v2, double gradient)
	{
		return v2 + (v1 - v2) * gradient;
	}

	static double gradient(PointF p0, PointF p1, PointF pn)
	{
		PointF p0p1 = p0 - p1;
		PointF pnp1 = pn - p1;

		if (abs(p0p1.x_) + abs(p0p1.y_) == 0)
			return 1;

		if (pnp1.x_ == p0p1.x_ && pnp1.y_ == p0p1.y_)
			return 1;

		return (abs(pnp1.x_) + abs(pnp1.y_)) / ( abs(p0p1.x_) + abs(p0p1.y_) );
	}
};


#endif
