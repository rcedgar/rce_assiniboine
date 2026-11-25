#pragma once

#include "naview.hpp"

// Write one figure: title, rotated naview drawing, wrapped dot-bracket,
// and the naview / RAG / RAG-dual graphs. Returns false if the file cannot
// be opened.
// wrap is the dot-bracket row length in symbols and must be >= 1.
// graph_width is the graph-row width as a fraction of the naview width and must be > 0.
// naview_only draws the structure and omits the title, dot-bracket, and graphs.
// The structure is still scaled to the dot-bracket row width.
// rainbow assigns one hue per group. Otherwise groups cycle through palette, or
// through the built-in 10-color palette when npalette is 0.
bool write_naview_svg(const char* path, const char* label, const char* sequence, const char* db,
                      bool circular, const NavLayout* layout, int wrap, double graph_width,
                      bool naview_only, bool rainbow, const char* const* palette, int npalette);

// True when text is a color SVG accepts: a keyword, #RGB[A], #RRGGBB[AA],
// rgb()/rgba(), or hsl()/hsla().
bool naview_svg_color(const char* text);
