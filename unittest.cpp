
#include "stdafx.h"

#include "unittest.h"

#include "engine3d/matrix.hpp"

Sample_1::Sample_1()
{
  // vertices for a cube
  pts_[0].set(1, 1, 1);
  pts_[1].set(-1, 1, 1);
  pts_[2].set(-1, -1, 1);
  pts_[3].set(1, -1, 1);
  pts_[4].set(1, 1, -1);
  pts_[5].set(-1, 1, -1);
  pts_[6].set(-1, -1, -1);
  pts_[7].set(1, -1, -1);

  rotate_degree_ = 0;
}

// paint cube vertices
void Sample_1::paintCubePoints(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, radian);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, radian);
  Utility::set_scale(scale_matrix, 20);
  Utility::set_translate(translate_matrix, 200, 200, 200);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  // get new position for each vertices
  for (int i = 0; i < 8; i++)
  {
    PointF pt = pts_[i] * world_matrix;

    // ignore the point if outside the screen
    if (pt.x_ <= 0 || pt.x_ >= width || pt.y_ <= 0 || pt.y_ >= height)
    {
      return;
    }

    // draw a white point RGB(255, 255, 255) on x and y
    int anchor = ((int)(pt.y_) * width + (int)(pt.x_)) * bytePerPixel;
    buf[anchor] = 0xFF;
    buf[anchor + 1] = 0xFF;
    buf[anchor + 2] = 0xFF;
  }

  return;
}

// paint cube with line
void Sample_1::paintCubeLines(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // 8 integer point
  PointI pts[8];

  // line points
  std::vector<PointI> lines;

  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, radian);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, radian);
  Utility::set_scale(scale_matrix, 20);
  Utility::set_translate(translate_matrix, 200, 200, 200);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  // get new position for each vertices
  for (int i = 0; i < 8; i++)
  {
    PointF pt = pts_[i] * world_matrix;
    pts[i] = pt;
  }

  // draw 12 lines
  Utility::bresenham3d(pts[0], pts[1], lines);
  Utility::bresenham3d(pts[1], pts[2], lines);
  Utility::bresenham3d(pts[2], pts[3], lines);
  Utility::bresenham3d(pts[0], pts[3], lines);

  Utility::bresenham3d(pts[4], pts[5], lines);
  Utility::bresenham3d(pts[5], pts[6], lines);
  Utility::bresenham3d(pts[6], pts[7], lines);
  Utility::bresenham3d(pts[4], pts[7], lines);

  Utility::bresenham3d(pts[0], pts[4], lines);
  Utility::bresenham3d(pts[1], pts[5], lines);
  Utility::bresenham3d(pts[2], pts[6], lines);
  Utility::bresenham3d(pts[3], pts[7], lines);

  for (PointI pti : lines)
  {
    // ignore the point if outside the screen
    if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
    {
      return;
    }

    // draw a white point RGB(255, 255, 255) on x and y
    int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
    buf[anchor] = 0xFF;
    buf[anchor + 1] = 0xFF;
    buf[anchor + 2] = 0xFF;
  }

  return;
}

Sample_2::Sample_2()
{
  // vertices for a face (triangle)
  face_.set(
    PointF(1, 0, 0),
    PointF(-1, 1, 0),
    PointF(-1, -1, 0)
  );

  rotate_degree_ = 0;
}

void Sample_2::paintFace(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // triagle points
  std::vector<PointF> points[3];

  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, radian);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, radian);
  Utility::set_scale(scale_matrix, 20);
  Utility::set_translate(translate_matrix, 200, 200, 200);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  // re-calculate the face vertices
  Face face = face_ * world_matrix;

  // calculate all points inside the face
  face.rasterize();

  // get new position for each points inside face
  bool drawed = false;
  for (size_t i = 0; i < face.points_.size(); i++)
  {
    PointI pti = face.points_.at(i);
    // skip the point that out of screen
    if (pti.x_ < 0 || pti.x_ > width || pti.y_ < 0 || pti.y_ > height)
    {
      continue;
    }

    // draw a white point RGB(255, 255, 255) on x and y
    int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
    buf[anchor] = 0xFF;
    buf[anchor + 1] = 0xFF;
    buf[anchor + 2] = 0xFF;

    drawed = true;
  }

  return;
}

Sample_3::Sample_3()
{
  // vertices for a cube
  mesh_.init("monkey.json");
  rotate_degree_ = 0;
}

void Sample_3::paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, -pi_/2.0);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, 0);
  Utility::set_scale(scale_matrix, 80);
  Utility::set_translate(translate_matrix, 500, 300, 200);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  for (Face fc : mesh_.faces_ )
  {
    // re-calculate the face vertices
    Face face = fc * world_matrix;

    // draw triangle's edges
    std::vector<PointI> lines;

    Utility::bresenham3d(face.pt1_, face.pt2_, lines);
    Utility::bresenham3d(face.pt1_, face.pt3_, lines);
    Utility::bresenham3d(face.pt2_, face.pt3_, lines);

    for (PointI pti : lines)
    {
      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      // draw a white point RGB(255, 255, 255) on x and y 
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      buf[anchor] = 0xFF;
      buf[anchor + 1] = 0xFF;
      buf[anchor + 2] = 0xFF;
    }

    lines.clear();
  }

  return;
}

Sample_4::Sample_4()
{
  // init 4 faces
  PointF pta(0, 0, 0);
  PointF ptb(1, 0, 0);
  PointF ptc(0, 1, 0);
  PointF ptd(0, 0, 1);

  abc_.set(pta, ptb, ptc);
  abd_.set(pta, ptb, ptd);
  acd_.set(pta, ptc, ptd);
  bcd_.set(ptb, ptc, ptd);

  n_abc_.set(0, 0, -1);
  n_abd_.set(0, -1, 0);
  n_acd_.set(-1, 0, 0);
  n_bcd_.set(1, 1, 1);

  rotate_degree_ = 0;
}

PointF Sample_4::center(Face face)
{
  return PointF(
      (face.pt1_.x_ + face.pt2_.x_ + face.pt3_.x_) / 3.0,
      (face.pt1_.y_ + face.pt2_.y_ + face.pt3_.y_) / 3.0,
      (face.pt1_.z_ + face.pt2_.z_ + face.pt3_.z_) / 3.0
    );
}

void Sample_4::paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // z buf
  int_fast32_t* zbuf = new int_fast32_t[width * height];
  
  memset(zbuf, 0, sizeof(int_fast32_t) * (width * height));

  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, radian);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, radian);
  Utility::set_scale(scale_matrix, 150);
  Utility::set_translate(translate_matrix, 500, 300, 600);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  // light position
  PointF light(500, 300, 1000);

  // 4 faces
  Face abc = abc_ * world_matrix;
  Face abd = abd_ * world_matrix;
  Face acd = acd_ * world_matrix;
  Face bcd = bcd_ * world_matrix;

  // the normal for 4 faces
  VectorF n_abc = VectorF::normalize(n_abc_ * world_matrix);
  VectorF n_abd = VectorF::normalize(n_abd_ * world_matrix);
  VectorF n_acd = VectorF::normalize(n_acd_ * world_matrix);
  VectorF n_bcd = VectorF::normalize(n_bcd_ * world_matrix);

  // vector from center of the faces and light
  VectorF abc_light = VectorF::normalize(light - center(abc));
  VectorF abd_light = VectorF::normalize(light - center(abd));
  VectorF acd_light = VectorF::normalize(light - center(acd));
  VectorF bcd_light = VectorF::normalize(light - center(bcd));

  // cos sita between light and normal
  double cos_abc_light = VectorF::cosine(abc_light, n_abc);
  double cos_abd_light = VectorF::cosine(abd_light, n_abd);
  double cos_acd_light = VectorF::cosine(acd_light, n_acd);
  double cos_bcd_light = VectorF::cosine(bcd_light, n_bcd);

  // draw abc
  if (cos_abc_light > 0)
  {
    abc.rasterize();

    for (PointI pti : abc.points_)
    {
      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      // ignore the point if there is other point in front of it
      if (pti.z_ < zbuf[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuf[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      buf[anchor] = (int)(255 * cos_abc_light);
      buf[anchor + 1] = (int)(255 * cos_abc_light);
      buf[anchor + 2] = (int)(255 * cos_abc_light);
    }
  }

  // draw abd
  if (cos_abd_light > 0)
  {
    abd.rasterize();

    for (PointI pti : abd.points_)
    {
      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      // ignore the point if there is other point in front of it
      if (pti.z_ < zbuf[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuf[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      buf[anchor] = (int)(255 * cos_abd_light);
      buf[anchor + 1] = (int)(255 * cos_abd_light);
      buf[anchor + 2] = (int)(255 * cos_abd_light);
    }
  }

  // draw acd
  if (cos_acd_light > 0)
  {
    acd.rasterize();

    for (PointI pti : acd.points_)
    {
      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      // ignore the point if there is other point in front of it
      if (pti.z_ < zbuf[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuf[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      buf[anchor] = (int)(255 * cos_acd_light);
      buf[anchor + 1] = (int)(255 * cos_acd_light);
      buf[anchor + 2] = (int)(255 * cos_acd_light);
    }
  }

  // draw bcd
  if (cos_bcd_light > 0)
  {
    bcd.rasterize();

    for (PointI pti : bcd.points_)
    {
      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      // ignore the point if there is other point in front of it
      if (pti.z_ < zbuf[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuf[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      buf[anchor] = (int)(255 * cos_bcd_light);
      buf[anchor + 1] = (int)(255 * cos_bcd_light);
      buf[anchor + 2] = (int)(255 * cos_bcd_light);
    }
  }
  else
  {
    int x = 0;
  }

  delete [] zbuf;

  return;
}

Sample_5::Sample_5()
{
  // vertices for a cube
  mesh_.init("monkey.json");
  rotate_degree_ = 0;
}

void Sample_5::paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // z buffer
  int_fast32_t* zbuffer = new int_fast32_t[width * height];

  // memset(zbuffer, 0, sizeof(int_fast32_t) * (width * height));
  for (size_t idx = 0; idx < (size_t)(width * height); idx++)
  {
    zbuffer[idx] = INT_FAST32_MIN;
  }

  // light 
  PointF light(900, 300, 1000);

  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, -pi_ / 2.0);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, 0);
  Utility::set_scale(scale_matrix, 200);
  Utility::set_translate(translate_matrix, 500, 300, 0);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  for (Face fc : mesh_.faces_)
  {
    // re-calculate the face vertices
    Face face = fc * world_matrix;
    std::vector<PointI> out;

    Utility::rasterize_triangle(
      light, 
      face.pt1_, face.pt1n_,
      face.pt2_, face.pt2n_,
      face.pt3_, face.pt3n_, 
      out);

    for ( size_t idx = 0; idx < out.size(); idx ++ )
    {
      PointI pti = out.at(idx);
      double cos = pti.cos_ * 1.0 / PointI::COS_PRECISENESS;

      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      if (pti.z_ < zbuffer[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuffer[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw a white point RGB(255, 255, 255) on x and y 
      int anchor = (pti.y_ * width + pti.x_) * bytePerPixel;

      if ( cos > 0 )
      {      
        buf[anchor] = (int)(255 * cos);
        buf[anchor + 1] = (int)(255 * cos);
        buf[anchor + 2] = (int)(255 * cos);
      }
    }
  }

  delete zbuffer;
  return;
}

Sample_6::Sample_6()
{
  // read bitmap
  bitmap_raw_ = NULL;
  std::ifstream fin("texture.bmp");

  if (fin.good())
  {
    unsigned char info[54];
    fin.read((char*)info, 54);

    bitmap_width_ = *(int32_t*)&info[18];
    bitmap_height_ = *(int32_t*)&info[22];

    int_fast32_t bitmap_size = 3 * bitmap_width_ * bitmap_height_;
    bitmap_raw_ = new unsigned char[bitmap_size];

    fin.read((char*)bitmap_raw_, bitmap_size);
  }

  // vertices for a cube
  PointF pta(0, 0, 0, 0, 0);
  PointF ptb(1, 0, 0, 399, 0);
  PointF ptc(0, 1, 0, 0, 399);
  PointF ptd(0, 0, 1, 399, 399);

  // non-accurate normal
  VectorF va(-1, -1, -1);
  VectorF vb(1, 0, 0);
  VectorF vc(0, 1, 0);
  VectorF vd(0, 0, 1);

  Face abc(pta, va, ptb, vb, ptc, vc );
  Face abd(pta, va, ptb, vd, ptd, vd );
  Face acd(pta, va, ptc, vc, ptd, vd );
  Face bcd(ptb, vb, ptc, vc, ptd, vd );

  mesh_.faces_.push_back(abc);
  mesh_.faces_.push_back(abd);
  mesh_.faces_.push_back(acd);
  mesh_.faces_.push_back(bcd);

  rotate_degree_ = 0;
}

Sample_6::~Sample_6()
{
  if (bitmap_raw_ != NULL)
  {
    delete[] bitmap_raw_;
  }
}

void Sample_6::paintMesh(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
{
  // z buffer
  int_fast32_t* zbuffer = new int_fast32_t[width * height];

  // memset(zbuffer, 0, sizeof(int_fast32_t) * (width * height));
  for (size_t idx = 0; idx < (size_t)(width * height); idx++)
  {
    zbuffer[idx] = INT_FAST32_MIN;
  }

  // light 
  PointF light(900, 300, 1000);

  // matrix
  Matrix<double> scale_matrix(4, 4);
  Matrix<double> translate_matrix(4, 4);
  Matrix<double> rotate_matrixX(4, 4);
  Matrix<double> rotate_matrixY(4, 4);
  Matrix<double> rotate_matrixZ(4, 4);

  // add one degree each time pain() has been called
  rotate_degree_ += 1;
  double radian = (rotate_degree_ % 360) * 2 * pi_ / 360;
  Utility::set_rotateX(rotate_matrixX, pi_ / 3.0);
  Utility::set_rotateY(rotate_matrixY, radian);
  Utility::set_rotateZ(rotate_matrixZ, 0);
  Utility::set_scale(scale_matrix, 200);
  Utility::set_translate(translate_matrix, 500, 300, 0);

  // final matrix
  Matrix<double> world_matrix =
    scale_matrix * rotate_matrixX * rotate_matrixY * rotate_matrixZ * translate_matrix;

  for (Face fc : mesh_.faces_)
  {
    // re-calculate the face vertices
    Face face = fc * world_matrix;
    std::vector<PointI> out;

    Utility::rasterize_triangle(
      light,
      face.pt1_, face.pt1n_,
      face.pt2_, face.pt2n_,
      face.pt3_, face.pt3n_,
      out);

    for (size_t idx = 0; idx < out.size(); idx++)
    {
      PointI pti = out.at(idx);
      double cos = pti.cos_ * 1.0 / PointI::COS_PRECISENESS;
      int_fast32_t u = pti.u_, v = pti.v_;

      // ignore the point if outside the screen 
      if (pti.x_ <= 0 || pti.x_ >= width || pti.y_ <= 0 || pti.y_ >= height)
      {
        continue;
      }

      if (pti.z_ < zbuffer[pti.y_ * width + pti.x_])
      {
        continue;
      }
      else
      {
        zbuffer[pti.y_ * width + pti.x_] = pti.z_;
      }

      // draw a white point RGB(255, 255, 255) on x and y 
      int buffer_anchor = (pti.y_ * width + pti.x_) * bytePerPixel;
      int bitmap_anchor = (v * bitmap_width_ + u) * 3;

      if (cos > 0)
      {
        if (bitmap_raw_ == NULL)
        {
          continue;
        }

        buf[buffer_anchor] = (int)(bitmap_raw_[bitmap_anchor] * cos);
        buf[buffer_anchor+1] = (int)(bitmap_raw_[bitmap_anchor+1] * cos);
        buf[buffer_anchor+2] = (int)(bitmap_raw_[bitmap_anchor+2] * cos);
      }
    }
  }

  delete zbuffer;
  return;
}