
#include "stdafx.h"

#include "unittest.h"

#include "engine3d/matrix.hpp"

Sample_1::Sample_1()
{
	// vertices for a cube
	pts_[0].set( 1,  1,  1);
	pts_[1].set(-1,  1,  1);
	pts_[2].set(-1, -1,  1);
	pts_[3].set( 1, -1,  1);
	pts_[4].set( 1,  1, -1);
	pts_[5].set(-1,  1, -1);
	pts_[6].set(-1, -1, -1);
	pts_[7].set( 1, -1, -1);

	rotate_degree_ = 0;
}

// paint
void Sample_1::paint(BYTE* buf, LONG width, LONG height, WORD bytePerPixel)
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
		buf[anchor+1] = 0xFF;
		buf[anchor+2] = 0xFF;
	}

	return;
}