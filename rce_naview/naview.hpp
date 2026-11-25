#pragma once

#include <string>
#include <vector>

// Secondary-structure layout from doc/naview_layout_algorithm.md.
// Coordinates use the reference frame: bond length 1, drawing root at the
// origin, +x right, +y up. Graph pictures use that same frame.

enum class NavKind {
    Exterior,
    Circle,
    Hairpin,
    Bulge,
    Interior,
    Multibranch,
    Helix5,
    Helix3
};

enum class NavSide { None, Five, Three };

struct NavPoint {
    double x;
    double y;
};

struct NavLine {
    double x0;
    double y0;
    double x1;
    double y1;
};

struct NavPicture {
    std::vector<NavPoint> vertices;
    std::vector<NavLine> lines;
};

struct NavLoop {
    NavKind kind;
    int start; // inclusive, or -1 when the closing loop has no unpaired base
    int end;
    int parent;         // -1 for loop 0
    int closing_helix;  // -1 for loop 0
    std::vector<int> child_helices;
};

struct NavHelix {
    int start5;
    int end5;
    int start3;
    int end3;
    int outer_loop;
    int inner_loop;
};

struct NavSegment {
    NavKind kind;
    int start;
    int end;
    int helix_id; // -1 when this is an unpaired segment
    int loop_id;  // -1 when this is a helix strand
    NavSide side;
};

struct NavLayout {
    std::vector<double> x;
    std::vector<double> y;
    std::vector<int> mate; // partner, or -1
    std::vector<NavLoop> loops;
    std::vector<NavHelix> helices;
    std::vector<NavSegment> segments;
    std::vector<int> segment_of; // per base
    std::string naview_graph;
    std::string rag_tree;
    std::string rag_dual;
    NavPicture naview_picture;
    NavPicture rag_picture;
    NavPicture dual_picture;
};

// db has length n and contains only '(', '.', ')'. Brackets must be balanced
// and properly nested. On rejection, *err (if non-null) points at a static
// message and the function returns false.
bool naview_layout(const char* db, int n, bool circular, NavLayout* out, const char** err);
