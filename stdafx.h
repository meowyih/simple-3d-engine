
#pragma once

#include "targetver.h"

#define WIN32_LEAN_AND_MEAN

// Windows header
#include <windows.h>

// C header
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <tchar.h>

// Other Header
#include "engine3d.h"
#include "engine3d/vectorf.hpp"
#include "engine3d/face.hpp"
#include "engine3d/pointf.hpp"
#include "engine3d/point.hpp"

VOID OnPaint(HDC hdc);
VOID face_test1(BYTE* buffer);