
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
        return;
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