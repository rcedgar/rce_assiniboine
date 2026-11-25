#include "svg_page.hpp"

#include "asserta.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

const double kMargin = 8.0;
const double kVGap = 10.0;
const double kTitlePx = 16.0;
const double kBondPx = 14.0;
const double kLetterPx = 9.0;
const double kPad = 8.0;
const double kDbPx = 12.0;
const double kAdvance = 7.2;
const double kRowPitch = 16.0;
const double kCaptionPx = 11.0;
const double kWash = 0.1;

struct Pt {
    double x;
    double y;
};

double cross(const Pt& o, const Pt& a, const Pt& b) {
    return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
}

std::vector<Pt> convex_hull(std::vector<Pt> pts) {
    std::sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b) {
        if (a.x < b.x) {
            return true;
        }
        if (a.x > b.x) {
            return false;
        }
        return a.y < b.y;
    });
    pts.erase(std::unique(pts.begin(), pts.end(),
                          [](const Pt& a, const Pt& b) { return a.x == b.x && a.y == b.y; }),
              pts.end());
    int n = static_cast<int>(pts.size());
    if (n <= 1) {
        return pts;
    }
    std::vector<Pt> h(static_cast<size_t>(2 * n));
    int k = 0;
    for (int i = 0; i < n; ++i) {
        while (k >= 2 && cross(h[static_cast<size_t>(k - 2)], h[static_cast<size_t>(k - 1)],
                               pts[static_cast<size_t>(i)]) <= 0.0) {
            k -= 1;
        }
        h[static_cast<size_t>(k++)] = pts[static_cast<size_t>(i)];
    }
    int lower = k;
    for (int i = n - 2; i >= 0; --i) {
        while (k > lower && cross(h[static_cast<size_t>(k - 2)], h[static_cast<size_t>(k - 1)],
                                  pts[static_cast<size_t>(i)]) <= 0.0) {
            k -= 1;
        }
        h[static_cast<size_t>(k++)] = pts[static_cast<size_t>(i)];
    }
    h.resize(static_cast<size_t>(k > 0 ? k - 1 : 0));
    return h;
}

void group_hsl(int index, int total, char* buf, size_t nbuf) {
    double hue = 0.0;
    if (total > 0) {
        hue = 360.0 * static_cast<double>(index) / static_cast<double>(total);
    }
    std::snprintf(buf, nbuf, "hsl(%.2f, 72%%, 40%%)", hue);
}

// Neighbors in this list are far apart in hue, and each is dark enough to read
// as a glyph on white. Groups past the tenth repeat from the start.
const char* kDefaultPalette[] = {
    "#c62828", "#1565c0", "#ef6c00", "#6a1b9a", "#2e7d32",
    "#ad1457", "#00838f", "#6d4c41", "#283593", "#9e9d24",
};
const int kDefaultPaletteN = 10;

void push_base(std::vector<int>* poly, int base) {
    if (poly->empty() || poly->back() != base) {
        poly->push_back(base);
    }
}

void emit_polygon(FILE* f, const std::vector<Pt>& pts, double ox, double oy, const char* color) {
    if (pts.size() < 3) {
        return;
    }
    std::fprintf(f, "<polygon points=\"");
    for (size_t i = 0; i < pts.size(); ++i) {
        std::fprintf(f, "%.2f,%.2f ", pts[i].x + ox, pts[i].y + oy);
    }
    std::fprintf(f, "\" fill=\"%s\" fill-opacity=\"%.2f\" stroke=\"none\"/>\n", color, kWash);
}

void emit_bases(FILE* f, const std::vector<Pt>& pts, const std::vector<int>& bases, double ox, double oy,
                const char* color) {
    if (bases.size() < 3) {
        return;
    }
    std::fprintf(f, "<polygon points=\"");
    for (size_t i = 0; i < bases.size(); ++i) {
        const Pt& p = pts[static_cast<size_t>(bases[i])];
        std::fprintf(f, "%.2f,%.2f ", p.x + ox, p.y + oy);
    }
    std::fprintf(f, "\" fill=\"%s\" fill-opacity=\"%.2f\" stroke=\"none\"/>\n", color, kWash);
}

// A band of total width 2*half along an open path. One base becomes a diamond.
void emit_band(FILE* f, const std::vector<Pt>& pts, const std::vector<int>& bases, double half, double ox,
               double oy, const char* color) {
    int m = static_cast<int>(bases.size());
    if (m < 1) {
        return;
    }
    std::vector<Pt> poly;
    if (m == 1) {
        const Pt& p = pts[static_cast<size_t>(bases[0])];
        Pt a;
        a.x = p.x - half;
        a.y = p.y;
        Pt b;
        b.x = p.x;
        b.y = p.y - half;
        Pt c;
        c.x = p.x + half;
        c.y = p.y;
        Pt d;
        d.x = p.x;
        d.y = p.y + half;
        poly.push_back(a);
        poly.push_back(b);
        poly.push_back(c);
        poly.push_back(d);
        emit_polygon(f, poly, ox, oy, color);
        return;
    }
    std::vector<Pt> left(static_cast<size_t>(m));
    std::vector<Pt> right(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) {
        const Pt& p = pts[static_cast<size_t>(bases[static_cast<size_t>(i)])];
        double dx = 0.0;
        double dy = 0.0;
        if (i + 1 < m) {
            const Pt& q = pts[static_cast<size_t>(bases[static_cast<size_t>(i + 1)])];
            dx += q.x - p.x;
            dy += q.y - p.y;
        }
        if (i > 0) {
            const Pt& q = pts[static_cast<size_t>(bases[static_cast<size_t>(i - 1)])];
            dx += p.x - q.x;
            dy += p.y - q.y;
        }
        double len = std::sqrt(dx * dx + dy * dy);
        if (len < 1e-9) {
            dx = 0.0;
            dy = 1.0;
            len = 1.0;
        }
        double nx = -dy / len;
        double ny = dx / len;
        left[static_cast<size_t>(i)].x = p.x + nx * half;
        left[static_cast<size_t>(i)].y = p.y + ny * half;
        right[static_cast<size_t>(i)].x = p.x - nx * half;
        right[static_cast<size_t>(i)].y = p.y - ny * half;
    }
    poly.insert(poly.end(), left.begin(), left.end());
    for (int i = m - 1; i >= 0; --i) {
        poly.push_back(right[static_cast<size_t>(i)]);
    }
    emit_polygon(f, poly, ox, oy, color);
}

void segment_bases(const NavSegment& seg, int n, std::vector<int>* bases) {
    bases->clear();
    if (seg.start <= seg.end) {
        for (int b = seg.start; b <= seg.end; ++b) {
            bases->push_back(b);
        }
        return;
    }
    for (int b = seg.start; b < n; ++b) {
        bases->push_back(b);
    }
    for (int b = 0; b <= seg.end; ++b) {
        bases->push_back(b);
    }
}

// Boundary of one loop. Child helices are crossed at their outer pair so the
// polygon does not cover the stem or anything nested inside it.
std::vector<int> loop_boundary(const NavLayout* layout, int loop_id, bool circular,
                               const std::vector<int>& start5_of) {
    const NavLoop& loop = layout->loops[static_cast<size_t>(loop_id)];
    int n = static_cast<int>(layout->x.size());
    std::vector<int> poly;
    auto child_end = [&](int base) -> int {
        if (base < 0 || base >= n) {
            return -1;
        }
        int hid = start5_of[static_cast<size_t>(base)];
        if (hid < 0) {
            return -1;
        }
        const NavHelix& h = layout->helices[static_cast<size_t>(hid)];
        if (h.outer_loop != loop_id) {
            return -1;
        }
        return h.end3;
    };
    if (loop.closing_helix >= 0) {
        const NavHelix& closing = layout->helices[static_cast<size_t>(loop.closing_helix)];
        int cur = closing.end5;
        int stop = closing.start3;
        push_base(&poly, cur);
        int guard = 0;
        while (cur != stop) {
            asserta(++guard <= n + 2);
            int nxt = cur + 1;
            asserta(nxt >= 0 && nxt < n);
            int jumped = child_end(nxt);
            if (jumped >= 0) {
                push_base(&poly, nxt);
                cur = jumped;
                push_base(&poly, cur);
            } else {
                cur = nxt;
                push_base(&poly, cur);
            }
        }
        return poly;
    }
    asserta(circular);
    int cur = 0;
    push_base(&poly, 0);
    int jumped0 = child_end(0);
    if (jumped0 >= 0) {
        cur = jumped0;
        push_base(&poly, cur);
    }
    int guard = 0;
    for (;;) {
        asserta(++guard <= n + 2);
        int nxt = cur + 1;
        if (nxt >= n) {
            break;
        }
        int jumped = child_end(nxt);
        if (jumped >= 0) {
            push_base(&poly, nxt);
            cur = jumped;
            push_base(&poly, cur);
        } else {
            cur = nxt;
            push_base(&poly, cur);
        }
    }
    return poly;
}

void append_escaped(std::string* out, const char* text, int n) {
    for (int i = 0; i < n; ++i) {
        unsigned char c = static_cast<unsigned char>(text[i]);
        if (c == '&') {
            out->append("&amp;");
        } else if (c == '<') {
            out->append("&lt;");
        } else if (c == '>') {
            out->append("&gt;");
        } else if (c == '"') {
            out->append("&quot;");
        } else if (c == '\'') {
            out->append("&apos;");
        } else {
            out->push_back(static_cast<char>(c));
        }
    }
}

void xml_text(FILE* f, const char* text) {
    if (!text) {
        return;
    }
    std::string esc;
    append_escaped(&esc, text, static_cast<int>(std::char_traits<char>::length(text)));
    std::fputs(esc.c_str(), f);
}

struct GraphBox {
    double width;
    double height;
    double minx;
    double miny;
    double maxx;
    double maxy;
    double scale;
};

GraphBox fit_graph(const NavPicture& pic, double target_h) {
    bool any = false;
    double minx = 0.0;
    double miny = 0.0;
    double maxx = 0.0;
    double maxy = 0.0;
    auto add = [&](double x, double y) {
        if (!any) {
            minx = maxx = x;
            miny = maxy = y;
            any = true;
            return;
        }
        if (x < minx) {
            minx = x;
        }
        if (x > maxx) {
            maxx = x;
        }
        if (y < miny) {
            miny = y;
        }
        if (y > maxy) {
            maxy = y;
        }
    };
    for (size_t i = 0; i < pic.vertices.size(); ++i) {
        add(pic.vertices[i].x, pic.vertices[i].y);
    }
    for (size_t i = 0; i < pic.lines.size(); ++i) {
        add(pic.lines[i].x0, pic.lines[i].y0);
        add(pic.lines[i].x1, pic.lines[i].y1);
    }
    if (!any) {
        minx = -0.5;
        maxx = 0.5;
        miny = -0.5;
        maxy = 0.5;
    }
    // Pad in graph units so a flat edge or a single point still has height.
    const double pad = 0.12;
    minx -= pad;
    maxx += pad;
    miny -= pad;
    maxy += pad;
    GraphBox box;
    box.minx = minx;
    box.miny = miny;
    box.maxx = maxx;
    box.maxy = maxy;
    double src_h = maxy - miny;
    double src_w = maxx - minx;
    asserta(src_h > 0.0 && src_w > 0.0);
    box.scale = target_h / src_h;
    box.height = target_h;
    box.width = src_w * box.scale;
    return box;
}

double page_x(const GraphBox& box, double left, double x) {
    return left + (x - box.minx) * box.scale;
}

double page_y(const GraphBox& box, double top, double y) {
    return top + (box.maxy - y) * box.scale;
}

}  // namespace

namespace {

const char* kSvgNames[] = {
    "aliceblue", "antiquewhite", "aqua", "aquamarine", "azure", "beige", "bisque", "black",
    "blanchedalmond", "blue", "blueviolet", "brown", "burlywood", "cadetblue", "chartreuse",
    "chocolate", "coral", "cornflowerblue", "cornsilk", "crimson", "currentcolor", "cyan",
    "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkgrey", "darkkhaki",
    "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred", "darksalmon",
    "darkseagreen", "darkslateblue", "darkslategray", "darkslategrey", "darkturquoise",
    "darkviolet", "deeppink", "deepskyblue", "dimgray", "dimgrey", "dodgerblue", "firebrick",
    "floralwhite", "forestgreen", "fuchsia", "gainsboro", "ghostwhite", "gold", "goldenrod",
    "gray", "green", "greenyellow", "grey", "honeydew", "hotpink", "indianred", "indigo",
    "ivory", "khaki", "lavender", "lavenderblush", "lawngreen", "lemonchiffon", "lightblue",
    "lightcoral", "lightcyan", "lightgoldenrodyellow", "lightgray", "lightgreen", "lightgrey",
    "lightpink", "lightsalmon", "lightseagreen", "lightskyblue", "lightslategray",
    "lightslategrey", "lightsteelblue", "lightyellow", "lime", "limegreen", "linen", "magenta",
    "maroon", "mediumaquamarine", "mediumblue", "mediumorchid", "mediumpurple", "mediumseagreen",
    "mediumslateblue", "mediumspringgreen", "mediumturquoise", "mediumvioletred", "midnightblue",
    "mintcream", "mistyrose", "moccasin", "navajowhite", "navy", "oldlace", "olive", "olivedrab",
    "orange", "orangered", "orchid", "palegoldenrod", "palegreen", "paleturquoise",
    "palevioletred", "papayawhip", "peachpuff", "peru", "pink", "plum", "powderblue", "purple",
    "rebeccapurple", "red", "rosybrown", "royalblue", "saddlebrown", "salmon", "sandybrown",
    "seagreen", "seashell", "sienna", "silver", "skyblue", "slateblue", "slategray", "slategrey",
    "snow", "springgreen", "steelblue", "tan", "teal", "thistle", "tomato", "transparent",
    "turquoise", "violet", "wheat", "white", "whitesmoke", "yellow", "yellowgreen",
};

bool hex_digit(char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
}

void skip_ws(const char** p) {
    while (**p == ' ' || **p == '\t') {
        *p += 1;
    }
}

bool take_number(const char** p, bool* percent) {
    skip_ws(p);
    const char* s = *p;
    if (*s == '+' || *s == '-') {
        s += 1;
    }
    bool digit = false;
    while (*s >= '0' && *s <= '9') {
        digit = true;
        s += 1;
    }
    if (*s == '.') {
        s += 1;
        while (*s >= '0' && *s <= '9') {
            digit = true;
            s += 1;
        }
    }
    if (!digit) {
        return false;
    }
    if (*s == 'e' || *s == 'E') {
        const char* exp = s + 1;
        if (*exp == '+' || *exp == '-') {
            exp += 1;
        }
        if (*exp < '0' || *exp > '9') {
            return false;
        }
        while (*exp >= '0' && *exp <= '9') {
            exp += 1;
        }
        s = exp;
    }
    *p = s;
    skip_ws(p);
    *percent = false;
    if (**p == '%') {
        *percent = true;
        *p += 1;
    }
    return true;
}

bool take_char(const char** p, char c) {
    skip_ws(p);
    if (**p != c) {
        return false;
    }
    *p += 1;
    return true;
}

bool fn_name(const char** p, const char* name) {
    const char* s = *p;
    while (*name) {
        char c = *s;
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
        if (c != *name) {
            return false;
        }
        s += 1;
        name += 1;
    }
    *p = s;
    return true;
}

// ncolor is 3. hue_first requires a plain hue and percent saturation and lightness.
// alpha, when set, is one extra component that may be a number or a percent.
bool color_fn(const char* p, const char* name, bool hue_first, bool alpha) {
    if (!fn_name(&p, name) || !take_char(&p, '(')) {
        return false;
    }
    bool first_percent = false;
    int n = 3 + (alpha ? 1 : 0);
    for (int i = 0; i < n; ++i) {
        if (i > 0 && !take_char(&p, ',')) {
            return false;
        }
        bool percent = false;
        if (!take_number(&p, &percent)) {
            return false;
        }
        if (alpha && i == 3) {
            continue;
        }
        if (hue_first && i == 0) {
            if (percent) {
                return false;
            }
        } else if (hue_first) {
            if (!percent) {
                return false;
            }
        } else if (i == 0) {
            first_percent = percent;
        } else if (percent != first_percent) {
            return false;
        }
    }
    if (!take_char(&p, ')')) {
        return false;
    }
    skip_ws(&p);
    return *p == '\0';
}

bool keyword_color(const char* text) {
    char lower[64];
    size_t n = 0;
    for (const char* p = text; *p; ++p) {
        if (n + 1 >= sizeof(lower)) {
            return false;
        }
        char c = *p;
        if (c >= 'A' && c <= 'Z') {
            c = static_cast<char>(c - 'A' + 'a');
        }
        lower[n++] = c;
    }
    lower[n] = '\0';
    for (size_t i = 0; i < sizeof(kSvgNames) / sizeof(kSvgNames[0]); ++i) {
        if (std::strcmp(lower, kSvgNames[i]) == 0) {
            return true;
        }
    }
    return false;
}

}  // namespace

bool naview_svg_color(const char* text) {
    if (!text) {
        return false;
    }
    while (*text == ' ' || *text == '\t') {
        text += 1;
    }
    const char* end = text;
    while (*end) {
        end += 1;
    }
    while (end > text && (end[-1] == ' ' || end[-1] == '\t')) {
        end -= 1;
    }
    if (end == text) {
        return false;
    }
    std::string body(text, end);
    if (body[0] == '#') {
        int nhex = 0;
        for (size_t i = 1; i < body.size(); ++i) {
            if (!hex_digit(body[i])) {
                return false;
            }
            nhex += 1;
        }
        return nhex == 3 || nhex == 4 || nhex == 6 || nhex == 8;
    }
    if (color_fn(body.c_str(), "rgba", false, true) || color_fn(body.c_str(), "rgb", false, false) ||
        color_fn(body.c_str(), "hsla", true, true) || color_fn(body.c_str(), "hsl", true, false)) {
        return true;
    }
    return keyword_color(body.c_str());
}

bool write_naview_svg(const char* path, const char* label, const char* sequence, const char* db,
                      bool circular, const NavLayout* layout, int wrap, double graph_width,
                      bool naview_only, bool rainbow, const char* const* palette, int npalette) {
    asserta(path && label && sequence && db && layout);
    asserta(wrap >= 1);
    asserta(graph_width > 0.0);
    asserta(npalette >= 0);
    asserta(npalette == 0 || palette);
    const int kWrap = wrap;
    int n = static_cast<int>(layout->x.size());
    asserta(n >= 1);
    asserta(static_cast<int>(layout->y.size()) == n);
    asserta(static_cast<int>(layout->segment_of.size()) == n);

    std::vector<Pt> pts(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        pts[static_cast<size_t>(i)].x = layout->x[static_cast<size_t>(i)];
        pts[static_cast<size_t>(i)].y = layout->y[static_cast<size_t>(i)];
    }
    std::vector<Pt> hull = convex_hull(pts);
    double best_h = 1.0e300;
    double best_w = -1.0;
    double best_ang = 0.0;
    bool have = false;
    int hn = static_cast<int>(hull.size());
    for (int i = 0; i < hn; ++i) {
        int j = (i + 1) % hn;
        double ang = std::atan2(hull[static_cast<size_t>(j)].y - hull[static_cast<size_t>(i)].y,
                                hull[static_cast<size_t>(j)].x - hull[static_cast<size_t>(i)].x);
        double co = std::cos(ang);
        double si = std::sin(ang);
        double minx = 0.0;
        double maxx = 0.0;
        double miny = 0.0;
        double maxy = 0.0;
        for (int k = 0; k < n; ++k) {
            double xr = pts[static_cast<size_t>(k)].x * co + pts[static_cast<size_t>(k)].y * si;
            double yr = -pts[static_cast<size_t>(k)].x * si + pts[static_cast<size_t>(k)].y * co;
            if (k == 0 || xr < minx) {
                minx = xr;
            }
            if (k == 0 || xr > maxx) {
                maxx = xr;
            }
            if (k == 0 || yr < miny) {
                miny = yr;
            }
            if (k == 0 || yr > maxy) {
                maxy = yr;
            }
        }
        double h = maxy - miny;
        double w = maxx - minx;
        bool better = !have || h < best_h - 1e-9 ||
                      (std::fabs(h - best_h) <= 1e-9 && w > best_w + 1e-9) ||
                      (std::fabs(h - best_h) <= 1e-9 && std::fabs(w - best_w) <= 1e-9 &&
                       std::fabs(ang) < std::fabs(best_ang) - 1e-12);
        if (better) {
            have = true;
            best_h = h;
            best_w = w;
            best_ang = ang;
        }
    }
    double co = std::cos(best_ang);
    double si = std::sin(best_ang);
    for (int i = 0; i < n; ++i) {
        double x = pts[static_cast<size_t>(i)].x;
        double y = pts[static_cast<size_t>(i)].y;
        pts[static_cast<size_t>(i)].x = (x * co + y * si) * kBondPx;
        pts[static_cast<size_t>(i)].y = (-x * si + y * co) * kBondPx;
    }
    double minx = pts[0].x;
    double maxx = pts[0].x;
    double miny = pts[0].y;
    double maxy = pts[0].y;
    for (int i = 1; i < n; ++i) {
        if (pts[static_cast<size_t>(i)].x < minx) {
            minx = pts[static_cast<size_t>(i)].x;
        }
        if (pts[static_cast<size_t>(i)].x > maxx) {
            maxx = pts[static_cast<size_t>(i)].x;
        }
        if (pts[static_cast<size_t>(i)].y < miny) {
            miny = pts[static_cast<size_t>(i)].y;
        }
        if (pts[static_cast<size_t>(i)].y > maxy) {
            maxy = pts[static_cast<size_t>(i)].y;
        }
    }
    for (int i = 0; i < n; ++i) {
        pts[static_cast<size_t>(i)].y = maxy - pts[static_cast<size_t>(i)].y;
    }
    double nav_w = (maxx - minx) + 2.0 * kPad;
    double nav_h = (maxy - miny) + 2.0 * kPad;
    for (int i = 0; i < n; ++i) {
        pts[static_cast<size_t>(i)].x += kPad - minx;
        pts[static_cast<size_t>(i)].y += kPad;
    }

    int row_chars = n < kWrap ? n : kWrap;
    double text_w = static_cast<double>(row_chars) * kAdvance;
    int nrows = (n + kWrap - 1) / kWrap;
    double text_h = static_cast<double>(nrows) * kRowPitch;

    // Fit to the dot-bracket row when the structure is wider. --naview-only
    // omits that row but keeps this scale, so the drawing matches the full figure.
    double nav_scale = 1.0;
    if (nav_w > text_w) {
        nav_scale = text_w / nav_w;
        for (int i = 0; i < n; ++i) {
            pts[static_cast<size_t>(i)].x *= nav_scale;
            pts[static_cast<size_t>(i)].y *= nav_scale;
        }
        nav_w *= nav_scale;
        nav_h *= nav_scale;
    }
    double letter_px = kLetterPx * nav_scale;
    double nav_stroke = 1.25 * nav_scale;
    double band_half = 0.35 * kBondPx * nav_scale;

    double graph_h = 0.0;
    double graph_gap = 0.0;
    GraphBox boxes[3];
    const NavPicture* pics[3] = {&layout->naview_picture, &layout->rag_picture, &layout->dual_picture};
    const char* captions[3] = {"naview", "RAG", "RAG dual"};
    const char* gids[3] = {"graph-naview", "graph-rag", "graph-rag-dual"};
    double graphs_w = 0.0;
    if (!naview_only) {
        graph_h = 0.25 * nav_h;
        if (graph_h < 1.0) {
            graph_h = 1.0;
        }
        graph_gap = 0.5 * graph_h;
        for (int i = 0; i < 3; ++i) {
            boxes[i] = fit_graph(*pics[i], graph_h);
            graphs_w += boxes[i].width;
        }
        graphs_w += 2.0 * graph_gap;
        asserta(graphs_w > 0.0);
        double gscale = (graph_width * nav_w) / graphs_w;
        for (int i = 0; i < 3; ++i) {
            boxes[i].scale *= gscale;
            boxes[i].width *= gscale;
            boxes[i].height *= gscale;
        }
        graph_gap *= gscale;
        graph_h *= gscale;
        graphs_w *= gscale;
    }

    // Dot-bracket at the top-left. Naview sits to its right, tops aligned.
    // Graphs are one row under the naview, as wide as graph_width * nav_w,
    // first graph flush with the naview's left edge.
    double text_left = kMargin;
    double nav_left = kMargin;
    double graphs_left = kMargin;
    double caption_h = 0.0;
    double y_title = kMargin + kTitlePx;
    double y_text = kMargin;
    double y_nav = kMargin;
    double y_graphs = kMargin;
    double bottom = kMargin + nav_h;
    double right = kMargin + nav_w;
    if (!naview_only) {
        caption_h = kCaptionPx + 4.0;
        nav_left = text_left + text_w + kVGap;
        graphs_left = nav_left;
        y_text = y_title + kVGap;
        y_nav = y_text;
        y_graphs = y_nav + nav_h + kVGap;
        double text_bottom = y_text + text_h;
        double graphs_bottom = y_graphs + graph_h + caption_h;
        bottom = text_bottom;
        if (graphs_bottom > bottom) {
            bottom = graphs_bottom;
        }
        right = nav_left + nav_w;
        if (graphs_left + graphs_w > right) {
            right = graphs_left + graphs_w;
        }
        int label_n = static_cast<int>(std::char_traits<char>::length(label));
        double title_w = 0.62 * kTitlePx * static_cast<double>(label_n);
        if (kMargin + title_w > right) {
            right = kMargin + title_w;
        }
    }
    bottom += kMargin;
    double page_w = right + kMargin;
    double page_h = bottom;

    //fprintf(stderr, "%s\n", path);
    FILE* f = std::fopen(path, "wb");
    if (!f) {
        return false;
    }
    std::fprintf(f,
                 "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
                 "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.2f\" height=\"%.2f\" "
                 "viewBox=\"0 0 %.2f %.2f\">\n",
                 page_w, page_h, page_w, page_h);

    if (!naview_only) {
        std::fprintf(f, "<g id=\"title\">\n");
        std::fprintf(f,
                     "<text x=\"%.2f\" y=\"%.2f\" font-family=\"sans-serif\" font-size=\"%.2f\" "
                     "fill=\"#111111\">",
                     kMargin, y_title, kTitlePx);
        xml_text(f, label);
        std::fprintf(f, "</text>\n</g>\n");
    }

    std::vector<int> helix_ci(layout->helices.size(), -1);
    std::vector<int> loop_ci(layout->loops.size(), -1);
    std::vector<int> base_ci(static_cast<size_t>(n), 0);
    int gcount = 0;
    for (int i = 0; i < n; ++i) {
        int sid = layout->segment_of[static_cast<size_t>(i)];
        asserta(sid >= 0 && sid < static_cast<int>(layout->segments.size()));
        const NavSegment& seg = layout->segments[static_cast<size_t>(sid)];
        int ci = -1;
        if (seg.helix_id >= 0) {
            asserta(seg.helix_id < static_cast<int>(layout->helices.size()));
            if (helix_ci[static_cast<size_t>(seg.helix_id)] < 0) {
                helix_ci[static_cast<size_t>(seg.helix_id)] = gcount;
                gcount += 1;
            }
            ci = helix_ci[static_cast<size_t>(seg.helix_id)];
        } else {
            asserta(seg.loop_id >= 0 && seg.loop_id < static_cast<int>(layout->loops.size()));
            if (loop_ci[static_cast<size_t>(seg.loop_id)] < 0) {
                loop_ci[static_cast<size_t>(seg.loop_id)] = gcount;
                gcount += 1;
            }
            ci = loop_ci[static_cast<size_t>(seg.loop_id)];
        }
        base_ci[static_cast<size_t>(i)] = ci;
    }
    asserta(gcount >= 1);
    std::vector<std::string> colors(static_cast<size_t>(gcount));
    for (int i = 0; i < gcount; ++i) {
        if (rainbow) {
            char buf[64];
            group_hsl(i, gcount, buf, sizeof(buf));
            colors[static_cast<size_t>(i)] = buf;
        } else if (npalette > 0) {
            colors[static_cast<size_t>(i)] = palette[i % npalette];
        } else {
            colors[static_cast<size_t>(i)] = kDefaultPalette[i % kDefaultPaletteN];
        }
    }
    std::vector<int> start5_of(static_cast<size_t>(n), -1);
    for (int h = 0; h < static_cast<int>(layout->helices.size()); ++h) {
        int s5 = layout->helices[static_cast<size_t>(h)].start5;
        asserta(s5 >= 0 && s5 < n);
        start5_of[static_cast<size_t>(s5)] = h;
    }

    std::fprintf(f, "<g id=\"naview\">\n");
    auto line = [&](double x0, double y0, double x1, double y1, const char* color, double width) {
        std::fprintf(f,
                     "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" "
                     "stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>\n",
                     x0, y0, x1, y1, color, width);
    };
    if (layout->helices.empty()) {
        asserta(!loop_ci.empty() && loop_ci[0] >= 0);
        const char* color = colors[static_cast<size_t>(loop_ci[0])].c_str();
        if (circular) {
            std::vector<int> all(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) {
                all[static_cast<size_t>(i)] = i;
            }
            if (n >= 3) {
                emit_bases(f, pts, all, nav_left, y_nav, color);
            } else {
                emit_band(f, pts, all, band_half, nav_left, y_nav, color);
            }
        } else {
            std::vector<int> bases;
            for (size_t s = 0; s < layout->segments.size(); ++s) {
                const NavSegment& seg = layout->segments[s];
                if (seg.loop_id != 0 || seg.helix_id >= 0) {
                    continue;
                }
                segment_bases(seg, n, &bases);
                emit_band(f, pts, bases, band_half, nav_left, y_nav, color);
            }
        }
    } else {
        for (int li = 0; li < static_cast<int>(layout->loops.size()); ++li) {
            if (loop_ci[static_cast<size_t>(li)] < 0) {
                continue;
            }
            const char* color = colors[static_cast<size_t>(loop_ci[static_cast<size_t>(li)])].c_str();
            if (li == 0 && !circular) {
                std::vector<int> bases;
                for (size_t s = 0; s < layout->segments.size(); ++s) {
                    const NavSegment& seg = layout->segments[s];
                    if (seg.loop_id != li || seg.helix_id >= 0) {
                        continue;
                    }
                    segment_bases(seg, n, &bases);
                    emit_band(f, pts, bases, band_half, nav_left, y_nav, color);
                }
            } else {
                std::vector<int> poly = loop_boundary(layout, li, circular, start5_of);
                emit_bases(f, pts, poly, nav_left, y_nav, color);
            }
        }
        for (int h = 0; h < static_cast<int>(layout->helices.size()); ++h) {
            const NavHelix& helix = layout->helices[static_cast<size_t>(h)];
            asserta(helix_ci[static_cast<size_t>(h)] >= 0);
            const char* color = colors[static_cast<size_t>(helix_ci[static_cast<size_t>(h)])].c_str();
            std::vector<int> poly;
            for (int k = helix.start5; k <= helix.end5; ++k) {
                poly.push_back(k);
            }
            for (int k = helix.end3; k >= helix.start3; --k) {
                poly.push_back(k);
            }
            emit_bases(f, pts, poly, nav_left, y_nav, color);
        }
    }
    for (int i = 0; i + 1 < n; ++i) {
        line(pts[static_cast<size_t>(i)].x + nav_left, pts[static_cast<size_t>(i)].y + y_nav,
             pts[static_cast<size_t>(i + 1)].x + nav_left, pts[static_cast<size_t>(i + 1)].y + y_nav,
             "#212121", nav_stroke);
    }
    if (circular && n > 1) {
        line(pts[static_cast<size_t>(n - 1)].x + nav_left, pts[static_cast<size_t>(n - 1)].y + y_nav,
             pts[0].x + nav_left, pts[0].y + y_nav, "#212121", nav_stroke);
    }
    for (int i = 0; i < n; ++i) {
        int j = layout->mate[static_cast<size_t>(i)];
        if (j > i) {
            const char* color = colors[static_cast<size_t>(base_ci[static_cast<size_t>(i)])].c_str();
            line(pts[static_cast<size_t>(i)].x + nav_left, pts[static_cast<size_t>(i)].y + y_nav,
                 pts[static_cast<size_t>(j)].x + nav_left, pts[static_cast<size_t>(j)].y + y_nav, color,
                 nav_stroke);
        }
    }
    for (int i = 0; i < n; ++i) {
        const char* color = colors[static_cast<size_t>(base_ci[static_cast<size_t>(i)])].c_str();
        std::fprintf(f,
                     "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" dominant-baseline=\"central\" "
                     "font-family=\"sans-serif\" font-size=\"%.2f\" fill=\"%s\">",
                     pts[static_cast<size_t>(i)].x + nav_left, pts[static_cast<size_t>(i)].y + y_nav,
                     letter_px, color);
        std::string one;
        append_escaped(&one, sequence + i, 1);
        std::fputs(one.c_str(), f);
        std::fprintf(f, "</text>\n");
    }
    std::fprintf(f, "</g>\n");

    if (!naview_only) {
    std::fprintf(f, "<g id=\"dotbracket\" font-family=\"ui-monospace, Consolas, monospace\" "
                    "font-size=\"%.2f\">\n",
                 kDbPx);
    for (int row = 0; row < nrows; ++row) {
        int begin = row * kWrap;
        int count = n - begin;
        if (count > kWrap) {
            count = kWrap;
        }
        double row_w = static_cast<double>(count) * kAdvance;
        double x0 = text_left + (text_w - row_w) * 0.5;
        double top = y_text + static_cast<double>(row) * kRowPitch;
        int at = 0;
        while (at < count) {
            int ci = base_ci[static_cast<size_t>(begin + at)];
            int end = at + 1;
            while (end < count && base_ci[static_cast<size_t>(begin + end)] == ci) {
                end += 1;
            }
            const char* color = colors[static_cast<size_t>(ci)].c_str();
            double x = x0 + static_cast<double>(at) * kAdvance;
            double w = static_cast<double>(end - at) * kAdvance;
            std::fprintf(f,
                         "<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\" "
                         "fill-opacity=\"%.2f\"/>\n",
                         x, top, w, kRowPitch, color, kWash);
            at = end;
        }
    }
    for (int row = 0; row < nrows; ++row) {
        int begin = row * kWrap;
        int count = n - begin;
        if (count > kWrap) {
            count = kWrap;
        }
        double row_w = static_cast<double>(count) * kAdvance;
        double x = text_left + (text_w - row_w) * 0.5;
        double y = y_text + static_cast<double>(row) * kRowPitch + kDbPx;
        std::fprintf(f, "<text x=\"%.2f\" y=\"%.2f\" textLength=\"%.2f\" lengthAdjust=\"spacingAndGlyphs\">",
                     x, y, row_w);
        int at = 0;
        while (at < count) {
            int ci = base_ci[static_cast<size_t>(begin + at)];
            int end = at + 1;
            while (end < count && base_ci[static_cast<size_t>(begin + end)] == ci) {
                end += 1;
            }
            std::fprintf(f, "<tspan fill=\"%s\">", colors[static_cast<size_t>(ci)].c_str());
            std::string run;
            append_escaped(&run, db + begin + at, end - at);
            std::fputs(run.c_str(), f);
            std::fprintf(f, "</tspan>");
            at = end;
        }
        std::fprintf(f, "</text>\n");
    }
    std::fprintf(f, "</g>\n");

    double gx = graphs_left;
    for (int g = 0; g < 3; ++g) {
        const GraphBox& box = boxes[g];
        double edge_px = box.scale;
        double gstroke = std::min(0.7, 0.15 * edge_px);
        double dot_r = std::min(1.1, 0.22 * edge_px);
        std::fprintf(f, "<g id=\"%s\">\n", gids[g]);
        const NavPicture& pic = *pics[g];
        for (size_t i = 0; i < pic.lines.size(); ++i) {
            const NavLine& L = pic.lines[i];
            line(page_x(box, gx, L.x0), page_y(box, y_graphs, L.y0), page_x(box, gx, L.x1),
                 page_y(box, y_graphs, L.y1), "#212121", gstroke);
        }
        for (size_t i = 0; i < pic.vertices.size(); ++i) {
            std::fprintf(f, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"#212121\"/>\n",
                         page_x(box, gx, pic.vertices[i].x), page_y(box, y_graphs, pic.vertices[i].y),
                         dot_r);
        }
        std::fprintf(f,
                     "<text x=\"%.2f\" y=\"%.2f\" text-anchor=\"middle\" font-family=\"sans-serif\" "
                     "font-size=\"%.2f\" fill=\"#212121\">%s</text>\n",
                     gx + box.width * 0.5, y_graphs + graph_h + kCaptionPx, kCaptionPx, captions[g]);
        std::fprintf(f, "</g>\n");
        gx += box.width + graph_gap;
    }
    }

    std::fprintf(f, "</svg>\n");
    int rc = std::fclose(f);
    return rc == 0;
}
