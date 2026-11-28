#pragma once

#include "naview.hpp"

// Write a 1-bit grayscale PNG of the structure, black strokes on white.
// width and height are pixels, each from 1 to 8192. The drawing uses the
// same min-height rotation as the SVG. Returns false if the file cannot
// be written.
bool write_naview_png(const char* path, const NavLayout* layout, bool circular, int width,
                      int height);
