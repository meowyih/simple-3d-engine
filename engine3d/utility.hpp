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
  inline static void bresenham2d(PointI pt1, PointI pt2, std::vector<PointI>& out)
  {
    int_fast32_t x0 = pt1.x_;
    int_fast32_t y0 = pt1.y_;
    int_fast32_t x1 = pt2.x_;
    int_fast32_t y1 = pt2.y_;
    int_fast32_t tmp;
    bool steep = abs(pt2.y_ - pt1.y_) > abs(pt2.x_ - pt1.x_) ? true : false;
    PointI from, to;
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
        PointI pt(y, x, pt1.z_);
        out.push_back(pt);
      }
      else
      {
        PointI pt(x, y, pt1.z_);
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

  inline static void bresenham2d2(PointI pt1, PointI pt2, std::vector<PointI>& out)
  {
    PointI pos(pt1);
    PointI delta = pt2 - pt1;
    PointI abs_delta(abs(delta.x_), abs(delta.y_), abs(delta.z_));
    PointI double_delta(abs_delta.x_ << 1, abs_delta.y_ << 1, abs_delta.z_ << 1);

    int_fast32_t x_inc, y_inc, z_inc;
    int_fast32_t cos_cur = pt1.cos_;
    int_fast32_t delta_cos = pt2.cos_ - pt1.cos_;

    x_inc = (delta.x_ < 0) ? -1 : 1;
    y_inc = (delta.y_ < 0) ? -1 : 1;
    z_inc = (delta.z_ < 0) ? -1 : 1;

    if (abs_delta.x_ >= abs_delta.y_ && abs_delta.x_ >= abs_delta.z_)
    {
      int_fast32_t err_1 = double_delta.y_ - abs_delta.x_;
      int_fast32_t err_2 = double_delta.z_ - abs_delta.x_;
      int_fast32_t cos_inc = ( abs_delta.x_ == 0 ) ? 0 : (delta_cos / abs_delta.x_);

      out.reserve(abs_delta.x_ + 1);

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
        cos_cur += cos_inc;
        pos.cos_ = cos_cur;
      }
    }
    else if (abs_delta.y_ >= abs_delta.x_ && abs_delta.y_ >= abs_delta.z_)
    {
      int_fast32_t err_1 = double_delta.x_ - abs_delta.y_;
      int_fast32_t err_2 = double_delta.z_ - abs_delta.y_;
      int_fast32_t cos_inc = (abs_delta.y_ == 0) ? 0 : (delta_cos / abs_delta.y_);

      out.reserve(abs_delta.y_ + 1);

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
        cos_cur += cos_inc;
        pos.cos_ = cos_cur;
      }
    }
    else
    {
      int_fast32_t err_1 = double_delta.y_ - abs_delta.z_;
      int_fast32_t err_2 = double_delta.x_ - abs_delta.z_;
      int_fast32_t cos_inc = (abs_delta.z_ == 0) ? 0 : (delta_cos / abs_delta.z_);

      // create a point to store previous PointI in out
      PointI lastpos(pt1);
      pt1.x_++; // just make is different from pt1

      (delta.x_ > delta.y_) ? out.reserve(abs_delta.x_ + 1) : out.reserve(abs_delta.y_ + 1);

      for (int_fast32_t i = 0; i < abs_delta.z_; i++)
      {
        if (lastpos.x_ == pt1.x_ && lastpos.y_ == pos.y_)
        {
          // if x, y are the same, we don't have to store the pixel
        }
        else
        {
          out.push_back(pos);
          lastpos = pos;
        }

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
        cos_cur += cos_inc;
        pos.cos_ = cos_cur;
      }
    }

    pos.cos_ = pt2.cos_;
    out.push_back(pos);
  }

  // line3d uses Bresenham's algorithm to generate the 3 dimensional points on a
  // line from (x1, y1, z1) to (x2, y2, z2), reference implementation found here
  // http://www.ict.griffith.edu.au/anthony/info/graphics/bresenham.procs (3D)
  inline static void bresenham3d(PointI pt1, PointI pt2, std::vector<PointI>& out)
  {
    PointI pos(pt1);
    PointI delta = pt2 - pt1;
    PointI abs_delta(abs(delta.x_), abs(delta.y_), abs(delta.z_));
    PointI double_delta(abs_delta.x_ << 1, abs_delta.y_ << 1, abs_delta.z_ << 1);

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

  inline static void rasterize_triangle(PointI pt1, PointI pt2, PointI pt3, std::vector<PointI>& out)
  {
    PointI top, middle, bottom;

    std::vector<PointI> vertices_top_middle;
    std::vector<PointI> vertices_top_bottom;
    std::vector<PointI> vertices_middle_bottom;

    // sort three points by x
    if (pt1.x_ >= pt2.x_ && pt1.x_ >= pt3.x_)
    {
      top = pt1;
      if (pt2.x_ >= pt3.x_)
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
    else if (pt2.x_ >= pt1.x_ && pt2.x_ >= pt3.x_)
    {
      top = pt2;
      if (pt1.x_ >= pt3.x_)
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
      if (pt1.x_ >= pt2.x_)
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

    // horizontal line 
    if (top.x_ == middle.x_ && middle.x_ == bottom.x_)
    {
      PointI pt_ymin, pt_ymax;

      if (top.y_ >= middle.y_ && top.y_ >= bottom.y_)
      {
        pt_ymax = top;
      }
      else if (middle.y_ >= bottom.y_)
      {
        pt_ymax = middle;
      }
      else
      {
        pt_ymax = bottom;
      }

      if (top.y_ <= middle.y_ && top.y_ <= bottom.y_)
      {
        pt_ymin = top;
      }
      else if (middle.y_ <= bottom.y_)
      {
        pt_ymin = middle;
      }
      else
      {
        pt_ymin = bottom;
      }

      bresenham2d2(pt_ymin, pt_ymax, out);

      return;
    }

    // calculate all points between three edges
    bresenham2d2(top, middle, vertices_top_middle);
    bresenham2d2(top, bottom, vertices_top_bottom);
    bresenham2d2(middle, bottom, vertices_middle_bottom);

    auto it1 = vertices_top_bottom.begin();
    auto it2 = vertices_top_middle.begin();
    auto it3 = vertices_middle_bottom.begin();

    int_fast32_t x_cur = it1->x_;

    // reserve space for result vector, this is only
    // for performance concern. Increase vector's size at 
    // runtime could impact performance.
    size_t reserve = abs(top.x_ - bottom.x_ + 1) * (max(abs(top.y_ - middle.y_), abs(top.y_ - bottom.y_)) + 1);

    if (reserve > 4096)
    {
      out.reserve(4096);
    }
    else
    {
      out.reserve(reserve);
    }

    // fill upper triangle if needed
    if (top.x_ != middle.x_)
    {
      while (it2 != vertices_top_middle.end() && it1 != vertices_top_bottom.end())
      {
        // calculate ymax and ymin in x_cur
        int_fast32_t ymax, ymin, zmax, zmin;
        if (it1->y_ > it2->y_)
        {
          ymax = it1->y_;
          ymin = it2->y_;
          zmax = it1->z_;
          zmin = it2->z_;
        }
        else
        {
          ymax = it2->y_;
          ymin = it1->y_;
          zmax = it2->z_;
          zmin = it1->z_;
        }

        // find next it1->x_ == it2->x_ == x_cur
        while (it1 != vertices_top_bottom.end())
        {
          if (it1->x_ == x_cur)
          {
            if (it1->y_ > ymax)
            {
              ymax = it1->y_;
              zmax = it1->z_;
            }
            if (it1->y_ < ymin)
            {
              ymin = it1->y_;
              zmin = it1->z_;
            }
          }
          else // if (it1->x_ == x_cur - 1)
          {
            break;
          }
          it1++;
        }

        while (it2 != vertices_top_middle.end())
        {
          if (it2->x_ == x_cur)
          {
            if (it2->y_ > ymax)
            {
              ymax = it2->y_;
              zmax = it2->z_;
            }
            if (it2->y_ < ymin)
            {
              ymin = it2->y_;
              zmin = it2->z_;
            }
          }
          else // if (it2->x_ == x_cur - 1)
          {
            break;
          }
          it2++;
        }

        std::vector<PointI> points;
        bresenham2d2(PointI(x_cur, ymin, zmin), PointI(x_cur, ymax, zmax), points);
        out.insert(out.end(), points.begin(), points.end());
        x_cur--;

      }
    }

    // fill lower triangle if needed
    if (middle.x_ != bottom.x_)
    {
      // the start point of it3 is the end of it2, which already drawed
      // let skip it
      while (it3 != vertices_middle_bottom.end() && it3->x_ != x_cur)
      {
        it3++;
      }

      while (it3 != vertices_middle_bottom.end() && it1 != vertices_top_bottom.end())
      {
        // calculate ymax and ymin in x_cur
        int_fast32_t ymax, ymin, zmax, zmin;
        if (it1->y_ > it3->y_)
        {
          ymax = it1->y_;
          ymin = it3->y_;
          zmax = it1->z_;
          zmin = it3->z_;
        }
        else
        {
          ymax = it3->y_;
          ymin = it1->y_;
          zmax = it3->z_;
          zmin = it1->z_;
        }

        // find next it1->x_ == it2->x_ == x_cur
        while (it1 != vertices_top_bottom.end())
        {
          if (it1->x_ == x_cur)
          {
            if (it1->y_ > ymax)
            {
              ymax = it1->y_;
              zmax = it1->z_;
            }
            if (it1->y_ < ymin)
            {
              ymin = it1->y_;
              zmin = it1->z_;
            }
          }
          else // if (it1->x_ == x_cur - 1)
          {
            break;
          }
          it1++;
        }

        while (it3 != vertices_middle_bottom.end())
        {
          if (it3->x_ == x_cur)
          {
            if (it3->y_ > ymax)
            {
              ymax = it3->y_;
              zmax = it3->z_;
            }
            if (it3->y_ < ymin)
            {
              ymin = it3->y_;
              zmin = it3->z_;
            }
          }
          else // if (it3->x_ == x_cur - 1)
          {
            break;
          }
          it3++;
        }

        std::vector<PointI> points;
        bresenham2d2(PointI(x_cur, ymin, zmin), PointI(x_cur, ymax, zmax), points);
        out.insert(out.end(), points.begin(), points.end());
        x_cur--;
      }
    }
  }

  inline static void rasterize_triangle(
    PointI light,
    PointI pt1, VectorF pt1n,
    PointI pt2, VectorF pt2n,
    PointI pt3, VectorF pt3n,
    std::vector<PointI>& out)
  {
    PointI top, middle, bottom;
    VectorF top_n, middle_n, bottom_n;
    double top_c, middle_c, bottom_c;

    std::vector<PointI> vertices_top_middle;
    std::vector<PointI> vertices_top_bottom;
    std::vector<PointI> vertices_middle_bottom;

    // sort three points by x
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

    // calculate the cos sita value between line and three vertices
    PointF delta_light_top(
      light.x_ - top.x_, 
      light.y_ - top.y_, 
      light.z_ - top.z_);
    PointF delta_light_middle(
      light.x_ - middle.x_, 
      light.y_ - middle.y_, 
      light.z_ - middle.z_);
    PointF delta_light_bottom(
      light.x_ - bottom.x_, 
      light.y_ - bottom.y_, 
      light.z_ - bottom.z_);
    VectorF top_light = VectorF::normalize(delta_light_top);
    VectorF middle_light = VectorF::normalize(delta_light_middle);
    VectorF bottom_light = VectorF::normalize(delta_light_bottom);

    top_c = VectorF::cosine(top_light, top_n);
    middle_c = VectorF::cosine(middle_light, middle_n);
    bottom_c = VectorF::cosine(bottom_light, bottom_n);

    top.cos_ = (int_fast32_t)(top_c * PointI::COS_PRECISENESS);
    middle.cos_ = (int_fast32_t)(middle_c * PointI::COS_PRECISENESS);
    bottom.cos_ = (int_fast32_t)(bottom_c * PointI::COS_PRECISENESS);

    // horizontal line
    if (top.x_ == middle.x_ && middle.x_ == bottom.x_)
    {
      PointI pt_ymin, pt_ymax;
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

      bresenham2d2(pt_ymin, pt_ymax, out);

      return;
    }

    // calculate all points between three edges
    bresenham2d2(static_cast<PointI>(top), static_cast<PointI>(middle), vertices_top_middle);
    bresenham2d2(static_cast<PointI>(top), static_cast<PointI>(bottom), vertices_top_bottom);
    bresenham2d2(static_cast<PointI>(middle), static_cast<PointI>(bottom), vertices_middle_bottom);

    auto it1 = vertices_top_bottom.begin();
    auto it2 = vertices_top_middle.begin();
    auto it3 = vertices_middle_bottom.begin();

    int_fast32_t x_cur = it1->x_;

    PointI topi(top), bottomi(bottom), middlei(middle);

    // reserve space for result vector, this is only
    // for performance concern. Increase vector's size at 
    // runtime could impact performance.
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

    // fill upper triangle, skip it if there is no upper triangle
    if (top.x_ != middle.x_)
    {
      while (it2 != vertices_top_middle.end() && it1 != vertices_top_bottom.end())
      {
        // calculate ymax and ymin in x_cur
        int_fast32_t ymax, ymin, zmax, zmin, cosmax, cosmin;
        bool ymax_is_it1;
        if (it1->y_ > it2->y_)
        {
          ymax = it1->y_;
          ymin = it2->y_;
          zmax = it1->z_;
          zmin = it2->z_;
          cosmax = it1->cos_;
          cosmin = it2->cos_;
        }
        else
        {
          ymax = it2->y_;
          ymin = it1->y_;
          zmax = it2->z_;
          zmin = it1->z_;
          cosmax = it2->cos_;
          cosmin = it1->cos_;
        }

        // find next it1->x_ == it2->x_ == x_cur
        while (it1 != vertices_top_bottom.end())
        {
          if (it1->x_ == x_cur)
          {
            if (it1->y_ > ymax)
            {
              ymax = it1->y_;
              zmax = it1->z_;
              cosmax = it1->cos_;
            }
            if (it1->y_ < ymin)
            {
              ymin = it1->y_;
              zmin = it1->z_;
              cosmin = it1->cos_;
            }
          }
          else // if (it1->x_ == x_cur - 1)
          {
            break;
          }
          it1++;
        }

        while (it2 != vertices_top_middle.end())
        {
          if (it2->x_ == x_cur)
          {
            if (it2->y_ > ymax)
            {
              ymax_is_it1 = false;
              ymax = it2->y_;
              zmax = it2->z_;
              cosmax = it2->cos_;
            }
            if (it2->y_ < ymin)
            {
              ymin = it2->y_;
              zmin = it2->z_;
              cosmin = it2->cos_;
            }
          }
          else // if (it2->x_ == x_cur - 1)
          {
            break;
          }
          it2++;
        }

        std::vector<PointI> points;
        bresenham2d2(PointI(x_cur, ymin, zmin, cosmin), PointI(x_cur, ymax, zmax, cosmax), points);
        out.insert(out.end(), points.begin(), points.end());
        x_cur--;
      }
    }

    // fill lower triangle, skip it if there is no upper triangle
    if (middle.x_ != bottom.x_)
    {
      // the start point of it3 is the end of it2, which already drawed
      // let skip it
      while (it3 != vertices_middle_bottom.end() && it3->x_ != x_cur)
      {
        it3++;
      }

      while (it3 != vertices_middle_bottom.end() && it1 != vertices_top_bottom.end())
      {
        // calculate ymax and ymin in x_cur
        int_fast32_t ymax, ymin, zmax, zmin, cosmax, cosmin;
        if (it1->y_ > it3->y_)
        {
          ymax = it1->y_;
          ymin = it3->y_;
          zmax = it1->z_;
          zmin = it3->z_;
          cosmax = it1->cos_;
          cosmin = it3->cos_;
        }
        else
        {
          ymax = it3->y_;
          ymin = it1->y_;
          zmax = it3->z_;
          zmin = it1->z_;
          cosmax = it3->cos_;
          cosmin = it1->cos_;
        }

        // find next it1->x_ == it2->x_ == x_cur
        while (it1 != vertices_top_bottom.end())
        {
          if (it1->x_ == x_cur)
          {
            if (it1->y_ > ymax)
            {
              ymax = it1->y_;
              zmax = it1->z_;
              cosmax = it1->cos_;
            }
            if (it1->y_ < ymin)
            {
              ymin = it1->y_;
              zmin = it1->z_;
              cosmin = it1->cos_;
            }
          }
          else // if (it1->x_ == x_cur - 1)
          {
            break;
          }
          it1++;
        }

        while (it3 != vertices_middle_bottom.end())
        {
          if (it3->x_ == x_cur)
          {
            if (it3->y_ > ymax)
            {
              ymax = it3->y_;
              zmax = it3->z_;
              cosmax = it3->cos_;
            }
            if (it3->y_ < ymin)
            {
              ymin = it3->y_;
              zmin = it3->z_;
              cosmin = it3->cos_;
            }
          }
          else // if (it3->x_ == x_cur - 1)
          {
            break;
          }
          it3++;
        }

        std::vector<PointI> points;
        bresenham2d2(PointI(x_cur, ymin, zmin, cosmin), PointI(x_cur, ymax, zmax, cosmax), points);
        out.insert(out.end(), points.begin(), points.end());
        x_cur--;
      }
    }

  }

  static double interpolation(double v1, double v2, double gradient)
  {
    return v2 + (v1 - v2) * gradient;
  }

   static double gradient(PointI p0, PointI p1, PointI pn)
  {
    PointI p0p1 = p0 - p1;
    PointI pnp1 = pn - p1;

    if (abs(p0p1.x_) + abs(p0p1.y_) == 0)
      return 1;

    if (pnp1.x_ == p0p1.x_ && pnp1.y_ == p0p1.y_)
      return 1;

    return (abs(pnp1.x_) + abs(pnp1.y_)) * 1.0 / (abs(p0p1.x_) + abs(p0p1.y_)) * 1.0;
  }
};



#endif
