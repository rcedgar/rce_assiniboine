#include "png_page.hpp"

#include "asserta.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

const double kBondPx = 14.0;
const int kMaxSide = 8192;

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

// Same rotation and y-down placement as write_naview_svg, before page padding.
void place_structure(const NavLayout* layout, std::vector<Pt>* pts) {
    int n = static_cast<int>(layout->x.size());
    pts->resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        (*pts)[static_cast<size_t>(i)].x = layout->x[static_cast<size_t>(i)];
        (*pts)[static_cast<size_t>(i)].y = layout->y[static_cast<size_t>(i)];
    }
    std::vector<Pt> hull = convex_hull(*pts);
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
            double xr = (*pts)[static_cast<size_t>(k)].x * co + (*pts)[static_cast<size_t>(k)].y * si;
            double yr = -(*pts)[static_cast<size_t>(k)].x * si + (*pts)[static_cast<size_t>(k)].y * co;
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
        double x = (*pts)[static_cast<size_t>(i)].x;
        double y = (*pts)[static_cast<size_t>(i)].y;
        (*pts)[static_cast<size_t>(i)].x = (x * co + y * si) * kBondPx;
        (*pts)[static_cast<size_t>(i)].y = (-x * si + y * co) * kBondPx;
    }
    double maxy = (*pts)[0].y;
    for (int i = 1; i < n; ++i) {
        if ((*pts)[static_cast<size_t>(i)].y > maxy) {
            maxy = (*pts)[static_cast<size_t>(i)].y;
        }
    }
    for (int i = 0; i < n; ++i) {
        (*pts)[static_cast<size_t>(i)].y = maxy - (*pts)[static_cast<size_t>(i)].y;
    }
}

struct Pix {
    int x;
    int y;
};

int clamp_px(long v, int limit) {
    if (v < 0) {
        return 0;
    }
    if (v >= limit) {
        return limit - 1;
    }
    return static_cast<int>(v);
}

void fit_pixels(const std::vector<Pt>& pts, int width, int height, std::vector<Pix>* out) {
    int n = static_cast<int>(pts.size());
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
    int margin = 2;
    if (width - 2 * margin < 1 || height - 2 * margin < 1) {
        margin = 0;
    }
    int inner_w = width - 2 * margin;
    int inner_h = height - 2 * margin;
    double xspan = inner_w > 1 ? static_cast<double>(inner_w - 1) : 0.0;
    double yspan = inner_h > 1 ? static_cast<double>(inner_h - 1) : 0.0;
    double bw = maxx - minx;
    double bh = maxy - miny;
    const double flat = 1.0e-9;
    double scale = 0.0;
    if (bw > flat && bh > flat) {
        double sx = xspan / bw;
        double sy = yspan / bh;
        scale = sx < sy ? sx : sy;
    } else if (bw > flat) {
        scale = xspan / bw;
    } else if (bh > flat) {
        scale = yspan / bh;
    }
    double drawn_w = bw * scale;
    double drawn_h = bh * scale;
    double ox = static_cast<double>(margin) + (xspan - drawn_w) * 0.5;
    double oy = static_cast<double>(margin) + (yspan - drawn_h) * 0.5;
    out->resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        long lx = std::lround(ox + (pts[static_cast<size_t>(i)].x - minx) * scale);
        long ly = std::lround(oy + (pts[static_cast<size_t>(i)].y - miny) * scale);
        (*out)[static_cast<size_t>(i)].x = clamp_px(lx, width);
        (*out)[static_cast<size_t>(i)].y = clamp_px(ly, height);
    }
}

void paint(std::vector<unsigned char>* bits, int width, int height, int stride, int x, int y) {
    if (x < 0 || y < 0 || x >= width || y >= height) {
        return;
    }
    unsigned char* row = bits->data() + static_cast<size_t>(y) * static_cast<size_t>(stride);
    unsigned char mask = static_cast<unsigned char>(1u << (7 - (x & 7)));
    row[x >> 3] = static_cast<unsigned char>(row[x >> 3] & static_cast<unsigned char>(~mask));
}

void stamp(std::vector<unsigned char>* bits, int width, int height, int stride, int x, int y) {
    for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
            paint(bits, width, height, stride, x + dx, y + dy);
        }
    }
}

void draw_line(std::vector<unsigned char>* bits, int width, int height, int stride, int x0, int y0,
               int x1, int y1) {
    int dx = std::abs(x1 - x0);
    int sx = x0 < x1 ? 1 : -1;
    int dy = -std::abs(y1 - y0);
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;
    while (true) {
        stamp(bits, width, height, stride, x0, y0);
        if (x0 == x1 && y0 == y1) {
            break;
        }
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
}

void clear_padding(std::vector<unsigned char>* bits, int width, int height, int stride) {
    int rem = width & 7;
    if (rem == 0) {
        return;
    }
    unsigned char mask = static_cast<unsigned char>(0xFFu << (8 - rem));
    for (int y = 0; y < height; ++y) {
        unsigned char* row = bits->data() + static_cast<size_t>(y) * static_cast<size_t>(stride);
        row[stride - 1] = static_cast<unsigned char>(row[stride - 1] & mask);
    }
}

void fill_crc_table(std::uint32_t table[256]) {
    for (std::uint32_t i = 0; i < 256; ++i) {
        std::uint32_t c = i;
        for (int k = 0; k < 8; ++k) {
            if ((c & 1u) != 0u) {
                c = 0xEDB88320u ^ (c >> 1);
            } else {
                c >>= 1;
            }
        }
        table[i] = c;
    }
}

std::uint32_t crc_bytes(const std::uint32_t table[256], std::uint32_t crc, const unsigned char* data,
                        size_t n) {
    for (size_t i = 0; i < n; ++i) {
        crc = table[(crc ^ data[i]) & 0xFFu] ^ (crc >> 8);
    }
    return crc;
}

void put_be32(unsigned char* dst, std::uint32_t v) {
    dst[0] = static_cast<unsigned char>((v >> 24) & 0xFFu);
    dst[1] = static_cast<unsigned char>((v >> 16) & 0xFFu);
    dst[2] = static_cast<unsigned char>((v >> 8) & 0xFFu);
    dst[3] = static_cast<unsigned char>(v & 0xFFu);
}

bool write_bytes(FILE* f, const void* data, size_t n) {
    if (n == 0) {
        return true;
    }
    return std::fwrite(data, 1, n, f) == n;
}

bool write_chunk(FILE* f, const std::uint32_t table[256], const char* type, const unsigned char* data,
                 size_t len) {
    unsigned char lenb[4];
    put_be32(lenb, static_cast<std::uint32_t>(len));
    if (!write_bytes(f, lenb, 4)) {
        return false;
    }
    if (!write_bytes(f, type, 4)) {
        return false;
    }
    if (!write_bytes(f, data, len)) {
        return false;
    }
    std::uint32_t crc = 0xFFFFFFFFu;
    crc = crc_bytes(table, crc, reinterpret_cast<const unsigned char*>(type), 4);
    if (len > 0) {
        crc = crc_bytes(table, crc, data, len);
    }
    crc ^= 0xFFFFFFFFu;
    unsigned char crcb[4];
    put_be32(crcb, crc);
    return write_bytes(f, crcb, 4);
}

std::uint32_t adler32(const unsigned char* data, size_t n) {
    std::uint32_t s1 = 1;
    std::uint32_t s2 = 0;
    for (size_t i = 0; i < n; ++i) {
        s1 = (s1 + data[i]) % 65521u;
        s2 = (s2 + s1) % 65521u;
    }
    return (s2 << 16) | s1;
}

// zlib wrapper around stored deflate blocks. No compression library.
void stored_zlib(const unsigned char* raw, size_t n, std::vector<unsigned char>* out) {
    out->clear();
    out->push_back(0x78);
    out->push_back(0x01);
    size_t off = 0;
    while (off < n) {
        size_t chunk = n - off;
        if (chunk > 65535u) {
            chunk = 65535u;
        }
        bool final = off + chunk == n;
        out->push_back(final ? 1 : 0);
        unsigned len = static_cast<unsigned>(chunk);
        unsigned nlen = (~len) & 0xFFFFu;
        out->push_back(static_cast<unsigned char>(len & 0xFFu));
        out->push_back(static_cast<unsigned char>((len >> 8) & 0xFFu));
        out->push_back(static_cast<unsigned char>(nlen & 0xFFu));
        out->push_back(static_cast<unsigned char>((nlen >> 8) & 0xFFu));
        out->insert(out->end(), raw + off, raw + off + chunk);
        off += chunk;
    }
    unsigned char sum[4];
    put_be32(sum, adler32(raw, n));
    out->insert(out->end(), sum, sum + 4);
}

bool write_png(FILE* f, const std::vector<unsigned char>& bits, int width, int height, int stride) {
    std::uint32_t table[256];
    fill_crc_table(table);
    const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    if (!write_bytes(f, sig, 8)) {
        return false;
    }
    unsigned char ihdr[13];
    put_be32(ihdr, static_cast<std::uint32_t>(width));
    put_be32(ihdr + 4, static_cast<std::uint32_t>(height));
    ihdr[8] = 1;
    ihdr[9] = 0;
    ihdr[10] = 0;
    ihdr[11] = 0;
    ihdr[12] = 0;
    if (!write_chunk(f, table, "IHDR", ihdr, 13)) {
        return false;
    }
    size_t row_bytes = static_cast<size_t>(stride) + 1;
    std::vector<unsigned char> raw(static_cast<size_t>(height) * row_bytes);
    for (int y = 0; y < height; ++y) {
        size_t dst = static_cast<size_t>(y) * row_bytes;
        raw[dst] = 0;
        std::memcpy(raw.data() + dst + 1,
                    bits.data() + static_cast<size_t>(y) * static_cast<size_t>(stride),
                    static_cast<size_t>(stride));
    }
    std::vector<unsigned char> zlib;
    stored_zlib(raw.data(), raw.size(), &zlib);
    if (!write_chunk(f, table, "IDAT", zlib.data(), zlib.size())) {
        return false;
    }
    return write_chunk(f, table, "IEND", 0, 0);
}

}  // namespace

bool write_naview_png(const char* path, const NavLayout* layout, bool circular, int width, int height) {
    asserta(path && layout);
    asserta(width >= 1 && width <= kMaxSide);
    asserta(height >= 1 && height <= kMaxSide);
    int n = static_cast<int>(layout->x.size());
    asserta(n >= 1);
    asserta(static_cast<int>(layout->y.size()) == n);
    asserta(static_cast<int>(layout->mate.size()) == n);

    std::vector<Pt> placed;
    place_structure(layout, &placed);
    std::vector<Pix> pix;
    fit_pixels(placed, width, height, &pix);

    int stride = (width + 7) / 8;
    std::vector<unsigned char> bits(static_cast<size_t>(stride) * static_cast<size_t>(height), 0xFF);
    for (int i = 0; i + 1 < n; ++i) {
        draw_line(&bits, width, height, stride, pix[static_cast<size_t>(i)].x,
                  pix[static_cast<size_t>(i)].y, pix[static_cast<size_t>(i + 1)].x,
                  pix[static_cast<size_t>(i + 1)].y);
    }
    if (circular && n > 1) {
        draw_line(&bits, width, height, stride, pix[static_cast<size_t>(n - 1)].x,
                  pix[static_cast<size_t>(n - 1)].y, pix[0].x, pix[0].y);
    }
    for (int i = 0; i < n; ++i) {
        int j = layout->mate[static_cast<size_t>(i)];
        if (j > i) {
            draw_line(&bits, width, height, stride, pix[static_cast<size_t>(i)].x,
                      pix[static_cast<size_t>(i)].y, pix[static_cast<size_t>(j)].x,
                      pix[static_cast<size_t>(j)].y);
        }
    }
    clear_padding(&bits, width, height, stride);

    FILE* f = std::fopen(path, "wb");
    if (!f) {
        return false;
    }
    bool ok = write_png(f, bits, width, height, stride);
    if (!ok) {
        int err = errno;
        std::fclose(f);
        std::remove(path);
        errno = err;
        return false;
    }
    if (std::fclose(f) != 0) {
        int err = errno;
        std::remove(path);
        errno = err;
        return false;
    }
    return true;
}
