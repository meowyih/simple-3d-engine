
#include "stdafx.h"

#include "unittest.h"

VOID  OnPaint(HDC hdc)
{
	// get screen/window information
	BITMAP structBitmapHeader;
	memset(&structBitmapHeader, 0, sizeof(BITMAP));

	HGDIOBJ hBitmap = GetCurrentObject(hdc, OBJ_BITMAP);
	GetObject(hBitmap, sizeof(BITMAP), &structBitmapHeader);

	// windows width and height
	LONG width = structBitmapHeader.bmWidth;
	LONG height = structBitmapHeader.bmHeight;

	// bits per pixel, normally it is 24 (3 bytes) or 32 (4 bytes)
	WORD bytesPerPixel = structBitmapHeader.bmBitsPixel / 8;

	//
	// allocate pixel buffer for drawing
	//
	BYTE* buffer = ScreenBuffer::get_instance().get(width, height, bytesPerPixel);
	LONG buffersize = width * height * bytesPerPixel;

	// set all pixels to black 
	for (int i = 0; i < buffersize; i = i + bytesPerPixel)
	{
		buffer[i] = 0;     // BLUE
		buffer[i + 1] = 0; // GREEN
		buffer[i + 2] = 0; // RED
	}

	// sample code 
	Sample_1::get_instance().paint( buffer, width, height, bytesPerPixel);

	// update the screen with buffer
	BITMAPINFOHEADER bih;
	bih.biSize = sizeof(BITMAPINFOHEADER);
	bih.biBitCount = structBitmapHeader.bmBitsPixel;
	bih.biClrImportant = 0;
	bih.biClrUsed = 0;
	bih.biCompression = BI_RGB;
	bih.biWidth = structBitmapHeader.bmWidth;
	bih.biHeight = structBitmapHeader.bmHeight;
	bih.biPlanes = 1;
	bih.biSizeImage = buffersize;

	BITMAPINFO bi;
	bi.bmiHeader = bih;
	SetDIBitsToDevice(
		hdc, 0, 0, 
		width, height, 
		0, 0, 0, height,
		buffer, &bi, DIB_RGB_COLORS);
}