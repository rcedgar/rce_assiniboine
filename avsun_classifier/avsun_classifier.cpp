// One-class avsun thumbnail classifier.
//
// The positive class is the largest single-linkage component at an identity-pose
// chamfer of 2 px. Other avsuns are discarded. That component is split by a
// tighter chamfer so near-copies stay in one fold. A low score means the drawing
// looks like a prototype from this class. Designs are ranked by validation
// recall at a 1% decoy false-positive rate.
//
// avsun_classifier --avsun DIR --decoy DIR --out DIR

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#include <windows.h>
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

#define asserta(cond)                                                                          \
    do {                                                                                       \
        if (!(cond)) {                                                                         \
            std::fprintf(stderr, "asserta failed: %s:%d: %s\n", __FILE__, __LINE__, #cond);   \
            std::abort();                                                                      \
        }                                                                                      \
    } while (0)

static void die(const char* message) {
    std::fprintf(stderr, "avsun_classifier: %s\n", message);
    std::exit(1);
}

static const float kClip = 8.f;
static const float kLink = 2.f;
static const int kOccW = 25;
static const int kOccH = 12;
static const int kOccN = kOccW * kOccH;
static const int kAng = 16;
static const int kRad = 8;
static const int kPolarN = kAng * kRad;
static const float kInf = 1.0e6f;

enum Method { kChamfer = 0, kCosine = 1, kPolar = 2 };

static const char* method_name(Method method) {
    if (method == kChamfer) {
        return "chamfer";
    }
    if (method == kCosine) {
        return "cosine";
    }
    return "polar";
}

struct Image {
    std::vector<uint16_t> xs;
    std::vector<uint16_t> ys;
    std::vector<float> dt;
    float occ[kOccN];
    float polar[kPolarN];
    int npts;
};

struct Pose {
    float dx;
    float dy;
    float rot;
    float scale;
};

struct Design {
    Method method;
    int k;
    int trans;
    int rot180;
    int scale;
    double val_recall_fpr1;
    int val_fp;
    double test_recall_fpr1;
    int test_fp;
    double test_recall_fpr5;
    double test_recall_fpr01;
    float threshold;
    int val_hits;
    int test_hits;
};

static int g_w = 0;
static int g_h = 0;
static std::vector<Image> g_im;

static void put_be32_read(const unsigned char* p, std::uint32_t* out) {
    *out = (std::uint32_t)p[0] << 24 | (std::uint32_t)p[1] << 16 | (std::uint32_t)p[2] << 8 |
           (std::uint32_t)p[3];
}

static std::uint32_t crc_update(const std::uint32_t table[256], std::uint32_t crc,
                                const unsigned char* data, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        crc = table[(crc ^ data[i]) & 0xffu] ^ (crc >> 8);
    }
    return crc;
}

static void fill_crc_table(std::uint32_t table[256]) {
    for (std::uint32_t i = 0; i < 256; ++i) {
        std::uint32_t c = i;
        for (int k = 0; k < 8; ++k) {
            if ((c & 1u) != 0u) {
                c = 0xedb88320u ^ (c >> 1);
            } else {
                c >>= 1;
            }
        }
        table[i] = c;
    }
}

static bool read_file(const char* path, std::vector<unsigned char>* out) {
    FILE* f = std::fopen(path, "rb");
    if (!f) {
        return false;
    }
    if (std::fseek(f, 0, SEEK_END) != 0) {
        std::fclose(f);
        return false;
    }
    long sz = std::ftell(f);
    if (sz < 0) {
        std::fclose(f);
        return false;
    }
    if (std::fseek(f, 0, SEEK_SET) != 0) {
        std::fclose(f);
        return false;
    }
    out->resize(static_cast<size_t>(sz));
    if (sz > 0 && std::fread(out->data(), 1, static_cast<size_t>(sz), f) != static_cast<size_t>(sz)) {
        std::fclose(f);
        return false;
    }
    std::fclose(f);
    return true;
}

// 1-bit grayscale PNG as written by rce_naview: stored zlib blocks, black is 0.
static bool decode_png(const std::vector<unsigned char>& file, std::vector<unsigned char>* ink, int* w,
                       int* h) {
    if (file.size() < 8) {
        return false;
    }
    const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    if (std::memcmp(file.data(), sig, 8) != 0) {
        return false;
    }
    std::uint32_t table[256];
    fill_crc_table(table);
    bool have_ihdr = false;
    int width = 0;
    int height = 0;
    std::vector<unsigned char> idat;
    size_t pos = 8;
    while (pos + 12 <= file.size()) {
        std::uint32_t len = 0;
        put_be32_read(file.data() + pos, &len);
        if (pos + 12 + static_cast<size_t>(len) > file.size()) {
            return false;
        }
        const unsigned char* type = file.data() + pos + 4;
        const unsigned char* data = file.data() + pos + 8;
        std::uint32_t crc = 0xffffffffu;
        crc = crc_update(table, crc, type, 4);
        if (len > 0) {
            crc = crc_update(table, crc, data, len);
        }
        crc ^= 0xffffffffu;
        std::uint32_t stored = 0;
        put_be32_read(file.data() + pos + 8 + len, &stored);
        if (crc != stored) {
            return false;
        }
        if (std::memcmp(type, "IHDR", 4) == 0) {
            if (len != 13) {
                return false;
            }
            std::uint32_t uw = 0;
            std::uint32_t uh = 0;
            put_be32_read(data, &uw);
            put_be32_read(data + 4, &uh);
            if (uw == 0 || uh == 0 || uw > 8192 || uh > 8192) {
                return false;
            }
            if (data[8] != 1 || data[9] != 0 || data[10] != 0 || data[11] != 0 || data[12] != 0) {
                return false;
            }
            width = static_cast<int>(uw);
            height = static_cast<int>(uh);
            have_ihdr = true;
        } else if (std::memcmp(type, "IDAT", 4) == 0) {
            idat.insert(idat.end(), data, data + len);
        } else if (std::memcmp(type, "IEND", 4) == 0) {
            break;
        }
        pos += 12 + static_cast<size_t>(len);
    }
    if (!have_ihdr || idat.size() < 6) {
        return false;
    }
    if (idat[0] != 0x78 || idat[1] != 0x01) {
        return false;
    }
    std::vector<unsigned char> raw;
    size_t ip = 2;
    bool final = false;
    while (!final) {
        if (ip + 5 > idat.size()) {
            return false;
        }
        unsigned char header = idat[ip];
        if (header != 0 && header != 1) {
            return false;
        }
        final = header == 1;
        unsigned len = (unsigned)idat[ip + 1] | ((unsigned)idat[ip + 2] << 8);
        unsigned nlen = (unsigned)idat[ip + 3] | ((unsigned)idat[ip + 4] << 8);
        if ((nlen ^ 0xffffu) != len) {
            return false;
        }
        ip += 5;
        if (ip + len > idat.size()) {
            return false;
        }
        raw.insert(raw.end(), idat.begin() + static_cast<std::ptrdiff_t>(ip),
                   idat.begin() + static_cast<std::ptrdiff_t>(ip + len));
        ip += len;
    }
    int stride = (width + 7) / 8;
    size_t need = static_cast<size_t>(height) * static_cast<size_t>(stride + 1);
    if (raw.size() < need) {
        return false;
    }
    ink->assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
    for (int y = 0; y < height; ++y) {
        const unsigned char* row = raw.data() + static_cast<size_t>(y) * static_cast<size_t>(stride + 1);
        if (row[0] != 0) {
            return false;
        }
        const unsigned char* bits = row + 1;
        for (int x = 0; x < width; ++x) {
            int bit = (bits[x >> 3] >> (7 - (x & 7))) & 1;
            if (bit == 0) {
                (*ink)[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)] = 1;
            }
        }
    }
    *w = width;
    *h = height;
    return true;
}

static void edt_1d(const float* f, float* d, int n, int* v, float* z) {
    int k = 0;
    v[0] = 0;
    z[0] = -1.0e20f;
    z[1] = 1.0e20f;
    for (int q = 1; q < n; ++q) {
        int vk = v[k];
        float s = ((f[q] + static_cast<float>(q * q)) - (f[vk] + static_cast<float>(vk * vk))) /
                  static_cast<float>(2 * q - 2 * vk);
        while (s <= z[k]) {
            k -= 1;
            asserta(k >= 0);
            vk = v[k];
            s = ((f[q] + static_cast<float>(q * q)) - (f[vk] + static_cast<float>(vk * vk))) /
                static_cast<float>(2 * q - 2 * vk);
        }
        k += 1;
        v[k] = q;
        z[k] = s;
        z[k + 1] = 1.0e20f;
    }
    k = 0;
    for (int q = 0; q < n; ++q) {
        while (z[k + 1] < static_cast<float>(q)) {
            k += 1;
        }
        int vk = v[k];
        float diff = static_cast<float>(q - vk);
        d[q] = diff * diff + f[vk];
    }
}

static void build_dt(const unsigned char* ink, int w, int h, std::vector<float>* dt) {
    const int n = w * h;
    std::vector<float> grid(static_cast<size_t>(n), kInf);
    for (int i = 0; i < n; ++i) {
        if (ink[i] != 0) {
            grid[static_cast<size_t>(i)] = 0.f;
        }
    }
    int side = w > h ? w : h;
    std::vector<float> f(static_cast<size_t>(side));
    std::vector<float> d(static_cast<size_t>(side));
    std::vector<int> v(static_cast<size_t>(side));
    std::vector<float> z(static_cast<size_t>(side) + 1);
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            f[static_cast<size_t>(y)] = grid[static_cast<size_t>(y * w + x)];
        }
        edt_1d(f.data(), d.data(), h, v.data(), z.data());
        for (int y = 0; y < h; ++y) {
            grid[static_cast<size_t>(y * w + x)] = d[static_cast<size_t>(y)];
        }
    }
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            f[static_cast<size_t>(x)] = grid[static_cast<size_t>(y * w + x)];
        }
        edt_1d(f.data(), d.data(), w, v.data(), z.data());
        for (int x = 0; x < w; ++x) {
            float dist = std::sqrt(d[static_cast<size_t>(x)]);
            if (dist > kClip) {
                dist = kClip;
            }
            grid[static_cast<size_t>(y * w + x)] = dist;
        }
    }
    *dt = std::move(grid);
}

static void self_test_edt() {
    const int w = 5;
    const int h = 5;
    unsigned char ink[25];
    std::memset(ink, 0, sizeof(ink));
    ink[2 * w + 2] = 1;
    std::vector<float> dt;
    build_dt(ink, w, h, &dt);
    asserta(dt[static_cast<size_t>(2 * w + 2)] == 0.f);
    float neighbor = dt[static_cast<size_t>(2 * w + 3)];
    float corner = dt[static_cast<size_t>(0)];
    asserta(std::fabs(neighbor - 1.f) < 1.0e-4f);
    asserta(std::fabs(corner - std::sqrt(8.f)) < 1.0e-3f);
}

static void fill_occ(const uint16_t* xs, const uint16_t* ys, int npts, int w, int h, float* occ) {
    std::memset(occ, 0, sizeof(float) * static_cast<size_t>(kOccN));
    for (int i = 0; i < npts; ++i) {
        int bx = static_cast<int>(xs[i]) * kOccW / w;
        int by = static_cast<int>(ys[i]) * kOccH / h;
        if (bx < 0) {
            bx = 0;
        }
        if (by < 0) {
            by = 0;
        }
        if (bx >= kOccW) {
            bx = kOccW - 1;
        }
        if (by >= kOccH) {
            by = kOccH - 1;
        }
        occ[by * kOccW + bx] += 1.f;
    }
}

static void fill_polar(const uint16_t* xs, const uint16_t* ys, int npts, float* hist) {
    std::memset(hist, 0, sizeof(float) * static_cast<size_t>(kPolarN));
    if (npts <= 0) {
        return;
    }
    double sx = 0.0;
    double sy = 0.0;
    for (int i = 0; i < npts; ++i) {
        sx += xs[i];
        sy += ys[i];
    }
    double cx = sx / npts;
    double cy = sy / npts;
    double rmax = 0.5 * std::hypot(static_cast<double>(g_w), static_cast<double>(g_h));
    if (rmax < 1.0) {
        rmax = 1.0;
    }
    for (int i = 0; i < npts; ++i) {
        double dx = static_cast<double>(xs[i]) - cx;
        double dy = static_cast<double>(ys[i]) - cy;
        double ang = std::atan2(dy, dx);
        if (ang < 0.0) {
            ang += 6.283185307179586;
        }
        int ia = static_cast<int>(ang * (static_cast<double>(kAng) / 6.283185307179586));
        if (ia < 0) {
            ia = 0;
        }
        if (ia >= kAng) {
            ia = kAng - 1;
        }
        double rad = std::hypot(dx, dy) / rmax;
        int ir = static_cast<int>(rad * kRad);
        if (ir < 0) {
            ir = 0;
        }
        if (ir >= kRad) {
            ir = kRad - 1;
        }
        hist[ia * kRad + ir] += 1.f;
    }
    float sum = 0.f;
    for (int i = 0; i < kPolarN; ++i) {
        sum += hist[i];
    }
    if (sum > 0.f) {
        float inv = 1.f / sum;
        for (int i = 0; i < kPolarN; ++i) {
            hist[i] *= inv;
        }
    }
}

static bool ends_with_png(const char* name) {
    size_t n = std::strlen(name);
    return n >= 4 && std::strcmp(name + n - 4, ".png") == 0;
}

static void list_pngs(const char* dir, std::vector<std::string>* out) {
    out->clear();
#if defined(_WIN32)
    std::string pattern = std::string(dir) + "\\*.png";
    WIN32_FIND_DATAA found;
    HANDLE handle = FindFirstFileA(pattern.c_str(), &found);
    if (handle == INVALID_HANDLE_VALUE) {
        die("cannot list png directory");
    }
    do {
        if ((found.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) == 0 && ends_with_png(found.cFileName)) {
            out->push_back(std::string(dir) + "\\" + found.cFileName);
        }
    } while (FindNextFileA(handle, &found));
    FindClose(handle);
#else
    DIR* dp = opendir(dir);
    if (!dp) {
        die("cannot list png directory");
    }
    while (dirent* ent = readdir(dp)) {
        if (ends_with_png(ent->d_name)) {
            out->push_back(std::string(dir) + "/" + ent->d_name);
        }
    }
    closedir(dp);
#endif
    std::sort(out->begin(), out->end());
    if (out->empty()) {
        die("png directory is empty");
    }
}

static void load_image(const char* path, Image* image) {
    std::vector<unsigned char> file;
    if (!read_file(path, &file)) {
        std::fprintf(stderr, "avsun_classifier: cannot read %s\n", path);
        std::exit(1);
    }
    std::vector<unsigned char> ink;
    int w = 0;
    int h = 0;
    if (!decode_png(file, &ink, &w, &h)) {
        std::fprintf(stderr, "avsun_classifier: bad png %s\n", path);
        std::exit(1);
    }
    if (g_w == 0) {
        g_w = w;
        g_h = h;
    }
    if (w != g_w || h != g_h) {
        die("thumbnail sizes differ");
    }
    image->xs.clear();
    image->ys.clear();
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            if (ink[static_cast<size_t>(y * w + x)] != 0) {
                image->xs.push_back(static_cast<uint16_t>(x));
                image->ys.push_back(static_cast<uint16_t>(y));
            }
        }
    }
    image->npts = static_cast<int>(image->xs.size());
    build_dt(ink.data(), w, h, &image->dt);
    fill_occ(image->xs.data(), image->ys.data(), image->npts, w, h, image->occ);
    fill_polar(image->xs.data(), image->ys.data(), image->npts, image->polar);
}

static void load_dir(const char* dir, int* begin) {
    std::vector<std::string> paths;
    list_pngs(dir, &paths);
    *begin = static_cast<int>(g_im.size());
    g_im.reserve(g_im.size() + paths.size());
    for (size_t i = 0; i < paths.size(); ++i) {
        Image image;
        load_image(paths[i].c_str(), &image);
        g_im.push_back(std::move(image));
        if ((i + 1) % 1000 == 0) {
            std::fprintf(stderr, "loaded %s %zu/%zu\n", dir, i + 1, paths.size());
        }
    }
    std::fprintf(stderr, "loaded %s %zu\n", dir, paths.size());
}

static float mean_id(const Image& src, const float* dt) {
    if (src.npts <= 0) {
        return kClip;
    }
    float sum = 0.f;
    const int w = g_w;
    for (int i = 0; i < src.npts; ++i) {
        sum += dt[static_cast<int>(src.ys[static_cast<size_t>(i)]) * w + src.xs[static_cast<size_t>(i)]];
    }
    return sum / static_cast<float>(src.npts);
}

static float mean_id_capped(const Image& src, const float* dt, float max_mean) {
    if (src.npts <= 0) {
        return max_mean + 1.f;
    }
    float limit = max_mean * static_cast<float>(src.npts);
    float sum = 0.f;
    const int w = g_w;
    for (int i = 0; i < src.npts; ++i) {
        sum += dt[static_cast<int>(src.ys[static_cast<size_t>(i)]) * w + src.xs[static_cast<size_t>(i)]];
        if (sum > limit) {
            return max_mean + 1.f;
        }
    }
    return sum / static_cast<float>(src.npts);
}

static float chamfer_id(int a, int b) {
    float ab = mean_id(g_im[static_cast<size_t>(a)], g_im[static_cast<size_t>(b)].dt.data());
    float ba = mean_id(g_im[static_cast<size_t>(b)], g_im[static_cast<size_t>(a)].dt.data());
    return 0.5f * (ab + ba);
}

static bool chamfer_le(int a, int b, float max_mean, float* dist) {
    float cap = max_mean * 2.f;
    float ab = mean_id_capped(g_im[static_cast<size_t>(a)], g_im[static_cast<size_t>(b)].dt.data(), cap);
    if (ab > cap) {
        return false;
    }
    float remain = cap - ab;
    float ba = mean_id_capped(g_im[static_cast<size_t>(b)], g_im[static_cast<size_t>(a)].dt.data(), remain);
    if (ba > remain) {
        return false;
    }
    *dist = 0.5f * (ab + ba);
    return true;
}

static float cosine_raw(const float* a, const float* b) {
    float dot = 0.f;
    float na = 0.f;
    float nb = 0.f;
    for (int i = 0; i < kOccN; ++i) {
        dot += a[i] * b[i];
        na += a[i] * a[i];
        nb += b[i] * b[i];
    }
    if (na <= 0.f || nb <= 0.f) {
        return 1.f;
    }
    return 1.f - dot / std::sqrt(na * nb);
}

static float polar_l1(const float* a, const float* b) {
    float best = 1.0e9f;
    for (int rs = -1; rs <= 1; ++rs) {
        for (int ash = -1; ash <= 1; ++ash) {
            float sum = 0.f;
            for (int ang = 0; ang < kAng; ++ang) {
                int src_ang = ang - ash;
                if (src_ang < 0) {
                    src_ang += kAng;
                }
                if (src_ang >= kAng) {
                    src_ang -= kAng;
                }
                for (int rad = 0; rad < kRad; ++rad) {
                    int src_rad = rad - rs;
                    float av = 0.f;
                    if (src_rad >= 0 && src_rad < kRad) {
                        av = a[src_ang * kRad + src_rad];
                    }
                    sum += std::fabs(av - b[ang * kRad + rad]);
                }
            }
            if (sum < best) {
                best = sum;
            }
        }
    }
    return best;
}

struct Work {
    std::vector<int> px;
    std::vector<int> py;
    std::vector<unsigned char> ink;
    std::vector<float> grid;
    std::vector<float> f;
    std::vector<float> d;
    std::vector<int> v;
    std::vector<float> z;
    float occ[kOccN];
    float polar[kPolarN];
};

static void prepare_work(Work* work) {
    int n = g_w * g_h;
    int side = g_w > g_h ? g_w : g_h;
    work->ink.assign(static_cast<size_t>(n), 0);
    work->grid.assign(static_cast<size_t>(n), 0.f);
    work->f.assign(static_cast<size_t>(side), 0.f);
    work->d.assign(static_cast<size_t>(side), 0.f);
    work->v.assign(static_cast<size_t>(side), 0);
    work->z.assign(static_cast<size_t>(side) + 1, 0.f);
}

static void pose_points(const Image& image, const Pose& pose, Work* work) {
    int n = image.npts;
    work->px.resize(static_cast<size_t>(n));
    work->py.resize(static_cast<size_t>(n));
    float cx = 0.5f * static_cast<float>(g_w - 1);
    float cy = 0.5f * static_cast<float>(g_h - 1);
    float rad = pose.rot * 0.01745329251f;
    float co = std::cos(rad);
    float si = std::sin(rad);
    for (int i = 0; i < n; ++i) {
        float x = static_cast<float>(image.xs[static_cast<size_t>(i)]) - cx;
        float y = static_cast<float>(image.ys[static_cast<size_t>(i)]) - cy;
        x *= pose.scale;
        y *= pose.scale;
        float xr = co * x - si * y + cx + pose.dx;
        float yr = si * x + co * y + cy + pose.dy;
        work->px[static_cast<size_t>(i)] = static_cast<int>(std::lround(xr));
        work->py[static_cast<size_t>(i)] = static_cast<int>(std::lround(yr));
    }
}

static void edt_buffer(Work* work) {
    const int w = g_w;
    const int h = g_h;
    const int n = w * h;
    for (int i = 0; i < n; ++i) {
        work->grid[static_cast<size_t>(i)] = work->ink[static_cast<size_t>(i)] != 0 ? 0.f : kInf;
    }
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            work->f[static_cast<size_t>(y)] = work->grid[static_cast<size_t>(y * w + x)];
        }
        edt_1d(work->f.data(), work->d.data(), h, work->v.data(), work->z.data());
        for (int y = 0; y < h; ++y) {
            work->grid[static_cast<size_t>(y * w + x)] = work->d[static_cast<size_t>(y)];
        }
    }
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            work->f[static_cast<size_t>(x)] = work->grid[static_cast<size_t>(y * w + x)];
        }
        edt_1d(work->f.data(), work->d.data(), w, work->v.data(), work->z.data());
        for (int x = 0; x < w; ++x) {
            float dist = std::sqrt(work->d[static_cast<size_t>(x)]);
            if (dist > kClip) {
                dist = kClip;
            }
            work->grid[static_cast<size_t>(y * w + x)] = dist;
        }
    }
}

static float mean_posed(const Work& work, int npts, const float* dt) {
    if (npts <= 0) {
        return kClip;
    }
    float sum = 0.f;
    for (int i = 0; i < npts; ++i) {
        int x = work.px[static_cast<size_t>(i)];
        int y = work.py[static_cast<size_t>(i)];
        if (static_cast<unsigned>(x) >= static_cast<unsigned>(g_w) ||
            static_cast<unsigned>(y) >= static_cast<unsigned>(g_h)) {
            sum += kClip;
        } else {
            sum += dt[y * g_w + x];
        }
    }
    return sum / static_cast<float>(npts);
}

static float mean_proto_on_query(const Image& proto, const float* query_dt) {
    if (proto.npts <= 0) {
        return kClip;
    }
    float sum = 0.f;
    for (int i = 0; i < proto.npts; ++i) {
        sum += query_dt[static_cast<int>(proto.ys[static_cast<size_t>(i)]) * g_w +
                        proto.xs[static_cast<size_t>(i)]];
    }
    return sum / static_cast<float>(proto.npts);
}

static void raster_posed(Work* work, int npts) {
    std::fill(work->ink.begin(), work->ink.end(), 0);
    for (int i = 0; i < npts; ++i) {
        int x = work->px[static_cast<size_t>(i)];
        int y = work->py[static_cast<size_t>(i)];
        if (static_cast<unsigned>(x) >= static_cast<unsigned>(g_w) ||
            static_cast<unsigned>(y) >= static_cast<unsigned>(g_h)) {
            continue;
        }
        work->ink[static_cast<size_t>(y * g_w + x)] = 1;
    }
    edt_buffer(work);
}

static void occ_posed(const Work& work, int npts, float* occ) {
    std::memset(occ, 0, sizeof(float) * static_cast<size_t>(kOccN));
    int inside = 0;
    for (int i = 0; i < npts; ++i) {
        int x = work.px[static_cast<size_t>(i)];
        int y = work.py[static_cast<size_t>(i)];
        if (static_cast<unsigned>(x) >= static_cast<unsigned>(g_w) ||
            static_cast<unsigned>(y) >= static_cast<unsigned>(g_h)) {
            continue;
        }
        int bx = x * kOccW / g_w;
        int by = y * kOccH / g_h;
        if (bx >= kOccW) {
            bx = kOccW - 1;
        }
        if (by >= kOccH) {
            by = kOccH - 1;
        }
        occ[by * kOccW + bx] += 1.f;
        inside += 1;
    }
    (void)inside;
}

static void polar_from_posed(const Work& work, int npts, float* hist) {
    std::vector<uint16_t> xs;
    std::vector<uint16_t> ys;
    xs.reserve(static_cast<size_t>(npts));
    ys.reserve(static_cast<size_t>(npts));
    for (int i = 0; i < npts; ++i) {
        int x = work.px[static_cast<size_t>(i)];
        int y = work.py[static_cast<size_t>(i)];
        if (static_cast<unsigned>(x) >= static_cast<unsigned>(g_w) ||
            static_cast<unsigned>(y) >= static_cast<unsigned>(g_h)) {
            continue;
        }
        xs.push_back(static_cast<uint16_t>(x));
        ys.push_back(static_cast<uint16_t>(y));
    }
    fill_polar(xs.data(), ys.data(), static_cast<int>(xs.size()), hist);
}

static float identity_dist(Method method, int a, int b) {
    if (method == kChamfer) {
        return chamfer_id(a, b);
    }
    if (method == kCosine) {
        return cosine_raw(g_im[static_cast<size_t>(a)].occ, g_im[static_cast<size_t>(b)].occ);
    }
    return polar_l1(g_im[static_cast<size_t>(a)].polar, g_im[static_cast<size_t>(b)].polar);
}

static float score_query(int query, const std::vector<int>& protos, const std::vector<Pose>& poses,
                         Method method, Work* work) {
    float best = 1.0e9f;
    const Image& image = g_im[static_cast<size_t>(query)];
    for (size_t p = 0; p < poses.size(); ++p) {
        const Pose& pose = poses[p];
        bool identity = pose.dx == 0.f && pose.dy == 0.f && pose.rot == 0.f && pose.scale == 1.f;
        if (method == kChamfer) {
            if (identity) {
                for (size_t i = 0; i < protos.size(); ++i) {
                    float dist = chamfer_id(query, protos[i]);
                    if (dist < best) {
                        best = dist;
                    }
                }
            } else {
                pose_points(image, pose, work);
                raster_posed(work, image.npts);
                for (size_t i = 0; i < protos.size(); ++i) {
                    const Image& proto = g_im[static_cast<size_t>(protos[i])];
                    float ab = mean_posed(*work, image.npts, proto.dt.data());
                    float ba = mean_proto_on_query(proto, work->grid.data());
                    float dist = 0.5f * (ab + ba);
                    if (dist < best) {
                        best = dist;
                    }
                }
            }
        } else if (method == kCosine) {
            const float* occ = image.occ;
            float posed_occ[kOccN];
            if (!identity) {
                pose_points(image, pose, work);
                occ_posed(*work, image.npts, posed_occ);
                occ = posed_occ;
            }
            for (size_t i = 0; i < protos.size(); ++i) {
                float dist = cosine_raw(occ, g_im[static_cast<size_t>(protos[i])].occ);
                if (dist < best) {
                    best = dist;
                }
            }
        } else {
            const float* hist = image.polar;
            float posed_hist[kPolarN];
            if (!identity) {
                pose_points(image, pose, work);
                polar_from_posed(*work, image.npts, posed_hist);
                hist = posed_hist;
            }
            for (size_t i = 0; i < protos.size(); ++i) {
                float dist = polar_l1(hist, g_im[static_cast<size_t>(protos[i])].polar);
                if (dist < best) {
                    best = dist;
                }
            }
        }
    }
    return best;
}

static int find_root(std::vector<int>* parent, int x) {
    int p = (*parent)[static_cast<size_t>(x)];
    if (p != x) {
        p = find_root(parent, p);
        (*parent)[static_cast<size_t>(x)] = p;
    }
    return p;
}

static void unite(std::vector<int>* parent, std::vector<int>* rank, int a, int b) {
    a = find_root(parent, a);
    b = find_root(parent, b);
    if (a == b) {
        return;
    }
    if ((*rank)[static_cast<size_t>(a)] < (*rank)[static_cast<size_t>(b)]) {
        (*parent)[static_cast<size_t>(a)] = b;
    } else if ((*rank)[static_cast<size_t>(a)] > (*rank)[static_cast<size_t>(b)]) {
        (*parent)[static_cast<size_t>(b)] = a;
    } else {
        (*parent)[static_cast<size_t>(b)] = a;
        (*rank)[static_cast<size_t>(a)] += 1;
    }
}

struct Link {
    int a;
    int b;
    float dist;
};

static int thread_count(int n) {
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) {
        threads = 1;
    }
    if (threads > 32) {
        threads = 32;
    }
    if (threads > n) {
        threads = n;
    }
    if (threads < 1) {
        threads = 1;
    }
    return threads;
}

static void collect_links(const std::vector<int>& ids, float max_mean, std::vector<Link>* links) {
    int n = static_cast<int>(ids.size());
    int threads = thread_count(n);
    std::vector<std::vector<Link>> parts(static_cast<size_t>(threads));
    std::atomic<int> seen(0);
    std::mutex print_mu;
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            for (int i = t; i < n; i += threads) {
                int ia = ids[static_cast<size_t>(i)];
                for (int j = i + 1; j < n; ++j) {
                    float dist = 0.f;
                    if (chamfer_le(ia, ids[static_cast<size_t>(j)], max_mean, &dist)) {
                        Link link;
                        link.a = i;
                        link.b = j;
                        link.dist = dist;
                        parts[static_cast<size_t>(t)].push_back(link);
                    }
                }
                int done = ++seen;
                if (done % 250 == 0) {
                    std::lock_guard<std::mutex> lock(print_mu);
                    std::fprintf(stderr, "link %.2f %d/%d\n", max_mean, done, n);
                }
            }
        });
    }
    for (size_t t = 0; t < pool.size(); ++t) {
        pool[t].join();
    }
    links->clear();
    for (int t = 0; t < threads; ++t) {
        const std::vector<Link>& part = parts[static_cast<size_t>(t)];
        links->insert(links->end(), part.begin(), part.end());
    }
}

static void groups_at(int n, const std::vector<Link>& links, float cut, std::vector<int>* group_of, int* n_groups,
                      int* largest) {
    std::vector<int> parent(static_cast<size_t>(n));
    std::vector<int> rank(static_cast<size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        parent[static_cast<size_t>(i)] = i;
    }
    for (size_t e = 0; e < links.size(); ++e) {
        if (links[e].dist <= cut) {
            unite(&parent, &rank, links[e].a, links[e].b);
        }
    }
    group_of->assign(static_cast<size_t>(n), 0);
    std::vector<int> root_id(static_cast<size_t>(n), -1);
    std::vector<int> sizes;
    int groups = 0;
    for (int i = 0; i < n; ++i) {
        int root = find_root(&parent, i);
        int id = root_id[static_cast<size_t>(root)];
        if (id < 0) {
            id = groups;
            root_id[static_cast<size_t>(root)] = id;
            sizes.push_back(0);
            groups += 1;
        }
        (*group_of)[static_cast<size_t>(i)] = id;
        sizes[static_cast<size_t>(id)] += 1;
    }
    int big = 0;
    for (int g = 0; g < groups; ++g) {
        if (sizes[static_cast<size_t>(g)] > big) {
            big = sizes[static_cast<size_t>(g)];
        }
    }
    *n_groups = groups;
    *largest = big;
}

static void largest_component(const std::vector<int>& group_of, std::vector<int>* kept, std::vector<int>* dropped) {
    int n = static_cast<int>(group_of.size());
    int groups = 0;
    for (int i = 0; i < n; ++i) {
        if (group_of[static_cast<size_t>(i)] + 1 > groups) {
            groups = group_of[static_cast<size_t>(i)] + 1;
        }
    }
    std::vector<int> sizes(static_cast<size_t>(groups), 0);
    for (int i = 0; i < n; ++i) {
        sizes[static_cast<size_t>(group_of[static_cast<size_t>(i)])] += 1;
    }
    int best = 0;
    for (int g = 1; g < groups; ++g) {
        if (sizes[static_cast<size_t>(g)] > sizes[static_cast<size_t>(best)]) {
            best = g;
        }
    }
    kept->clear();
    dropped->clear();
    for (int i = 0; i < n; ++i) {
        if (group_of[static_cast<size_t>(i)] == best) {
            kept->push_back(i);
        } else {
            dropped->push_back(i);
        }
    }
}

struct Rng {
    std::uint32_t state;
    std::uint32_t next() {
        std::uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return x;
    }
};

static bool split_groups(const std::vector<int>& ids, const std::vector<int>& group_of, std::vector<int>* train,
                         std::vector<int>* val, std::vector<int>* test) {
    int n = static_cast<int>(ids.size());
    asserta(static_cast<int>(group_of.size()) == n);
    train->clear();
    val->clear();
    test->clear();
    int groups = 0;
    for (int i = 0; i < n; ++i) {
        if (group_of[static_cast<size_t>(i)] + 1 > groups) {
            groups = group_of[static_cast<size_t>(i)] + 1;
        }
    }
    std::vector<std::vector<int>> members(static_cast<size_t>(groups));
    for (int i = 0; i < n; ++i) {
        members[static_cast<size_t>(group_of[static_cast<size_t>(i)])].push_back(ids[static_cast<size_t>(i)]);
    }
    std::vector<int> order(static_cast<size_t>(groups));
    for (int g = 0; g < groups; ++g) {
        order[static_cast<size_t>(g)] = g;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        int sa = static_cast<int>(members[static_cast<size_t>(a)].size());
        int sb = static_cast<int>(members[static_cast<size_t>(b)].size());
        if (sa != sb) {
            return sa > sb;
        }
        return members[static_cast<size_t>(a)][0] < members[static_cast<size_t>(b)][0];
    });
    int target[3] = {(n + 1) / 2, (n + 2) / 4, n - ((n + 1) / 2) - ((n + 2) / 4)};
    int count[3] = {0, 0, 0};
    std::vector<int>* dest[3] = {train, val, test};
    for (int oi = 0; oi < groups; ++oi) {
        int g = order[static_cast<size_t>(oi)];
        int best = 0;
        int best_room = target[0] - count[0];
        for (int b = 1; b < 3; ++b) {
            int room = target[b] - count[b];
            if (room > best_room) {
                best = b;
                best_room = room;
            }
        }
        const std::vector<int>& mem = members[static_cast<size_t>(g)];
        dest[best]->insert(dest[best]->end(), mem.begin(), mem.end());
        count[best] += static_cast<int>(mem.size());
    }
    return !train->empty() && !val->empty() && !test->empty();
}

static void split_decoys(int begin, int n, std::vector<int>* val, std::vector<int>* test) {
    std::vector<int> ids(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        ids[static_cast<size_t>(i)] = begin + i;
    }
    Rng rng;
    rng.state = 1u;
    for (int i = n - 1; i > 0; --i) {
        int j = static_cast<int>(rng.next() % static_cast<std::uint32_t>(i + 1));
        std::swap(ids[static_cast<size_t>(i)], ids[static_cast<size_t>(j)]);
    }
    int nval = n / 2;
    val->assign(ids.begin(), ids.begin() + nval);
    test->assign(ids.begin() + nval, ids.end());
    std::fprintf(stderr, "split decoy val %zu test %zu\n", val->size(), test->size());
}

static int pick_seed(Method method, const std::vector<int>& train) {
    int n = static_cast<int>(train.size());
    asserta(n > 0);
    int refs = 128;
    if (refs > n) {
        refs = n;
    }
    int best = 0;
    double best_sum = 1.0e300;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int r = 0; r < refs; ++r) {
            int j = (r * n) / refs;
            if (j == i) {
                continue;
            }
            sum += identity_dist(method, train[static_cast<size_t>(i)], train[static_cast<size_t>(j)]);
        }
        if (sum < best_sum) {
            best_sum = sum;
            best = i;
        }
    }
    return best;
}

static float medoid_radius(const std::vector<int>& ids) {
    int at = pick_seed(kChamfer, ids);
    int medoid = ids[static_cast<size_t>(at)];
    float radius = 0.f;
    for (size_t i = 0; i < ids.size(); ++i) {
        float dist = chamfer_id(medoid, ids[i]);
        if (dist > radius) {
            radius = dist;
        }
    }
    return radius;
}

static std::vector<int> farthest_prototypes(Method method, const std::vector<int>& train, int k) {
    int n = static_cast<int>(train.size());
    if (k > n) {
        k = n;
    }
    asserta(k >= 1);
    int seed_at = pick_seed(method, train);
    std::vector<int> proto;
    proto.push_back(train[static_cast<size_t>(seed_at)]);
    std::vector<float> nearest(static_cast<size_t>(n), 1.0e9f);
    std::vector<char> used(static_cast<size_t>(n), 0);
    used[static_cast<size_t>(seed_at)] = 1;
    nearest[static_cast<size_t>(seed_at)] = 0.f;
    for (int t = 1; t < k; ++t) {
        int last = proto.back();
        for (int i = 0; i < n; ++i) {
            if (used[static_cast<size_t>(i)] != 0) {
                continue;
            }
            float dist = identity_dist(method, train[static_cast<size_t>(i)], last);
            if (dist < nearest[static_cast<size_t>(i)]) {
                nearest[static_cast<size_t>(i)] = dist;
            }
        }
        int best_i = -1;
        float best_d = -1.f;
        for (int i = 0; i < n; ++i) {
            if (used[static_cast<size_t>(i)] != 0) {
                continue;
            }
            if (nearest[static_cast<size_t>(i)] > best_d) {
                best_d = nearest[static_cast<size_t>(i)];
                best_i = i;
            }
        }
        asserta(best_i >= 0);
        used[static_cast<size_t>(best_i)] = 1;
        proto.push_back(train[static_cast<size_t>(best_i)]);
    }
    return proto;
}

static void make_poses(int trans, int rot180, int scale, std::vector<Pose>* poses) {
    poses->clear();
    float shifts[3] = {0.f, -4.f, 4.f};
    int nshift = trans ? 3 : 1;
    float rots[2] = {0.f, 180.f};
    int nrot = rot180 ? 2 : 1;
    float scales[3] = {1.f, 0.92f, 1.08f};
    int nscale = scale ? 3 : 1;
    for (int is = 0; is < nscale; ++is) {
        for (int ir = 0; ir < nrot; ++ir) {
            for (int iy = 0; iy < nshift; ++iy) {
                for (int ix = 0; ix < nshift; ++ix) {
                    Pose pose;
                    pose.dx = shifts[ix];
                    pose.dy = shifts[iy];
                    pose.rot = rots[ir];
                    pose.scale = scales[is];
                    poses->push_back(pose);
                }
            }
        }
    }
}

static void score_list(const std::vector<int>& ids, const std::vector<int>& protos,
                       const std::vector<Pose>& poses, Method method, std::vector<float>* scores) {
    scores->assign(ids.size(), 0.f);
    int threads = static_cast<int>(std::thread::hardware_concurrency());
    if (threads < 1) {
        threads = 1;
    }
    if (threads > 32) {
        threads = 32;
    }
    if (threads > static_cast<int>(ids.size())) {
        threads = static_cast<int>(ids.size());
    }
    if (threads < 1) {
        return;
    }
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            Work work;
            prepare_work(&work);
            for (size_t i = static_cast<size_t>(t); i < ids.size(); i += static_cast<size_t>(threads)) {
                (*scores)[i] = score_query(ids[i], protos, poses, method, &work);
            }
        });
    }
    for (size_t t = 0; t < pool.size(); ++t) {
        pool[t].join();
    }
}

static float threshold_at_fpr(std::vector<float> scores, double fpr) {
    if (scores.empty()) {
        return 0.f;
    }
    std::sort(scores.begin(), scores.end());
    int n = static_cast<int>(scores.size());
    int k = static_cast<int>(std::floor(fpr * static_cast<double>(n) + 1.0e-12));
    if (k <= 0) {
        return scores[0] - 1.f;
    }
    if (k > n) {
        k = n;
    }
    return scores[static_cast<size_t>(k - 1)];
}

static double recall_le(const std::vector<float>& scores, float threshold) {
    if (scores.empty()) {
        return 0.0;
    }
    int hits = 0;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] <= threshold) {
            hits += 1;
        }
    }
    return static_cast<double>(hits) / static_cast<double>(scores.size());
}

static int count_le(const std::vector<float>& scores, float threshold) {
    int hits = 0;
    for (size_t i = 0; i < scores.size(); ++i) {
        if (scores[i] <= threshold) {
            hits += 1;
        }
    }
    return hits;
}

static void list_components(const std::vector<int>& group_of, std::vector<std::vector<int>>* comps) {
    int n = static_cast<int>(group_of.size());
    int groups = 0;
    for (int i = 0; i < n; ++i) {
        if (group_of[static_cast<size_t>(i)] + 1 > groups) {
            groups = group_of[static_cast<size_t>(i)] + 1;
        }
    }
    comps->clear();
    comps->resize(static_cast<size_t>(groups));
    for (int i = 0; i < n; ++i) {
        (*comps)[static_cast<size_t>(group_of[static_cast<size_t>(i)])].push_back(i);
    }
    std::sort(comps->begin(), comps->end(), [](const std::vector<int>& a, const std::vector<int>& b) {
        if (a.size() != b.size()) {
            return a.size() > b.size();
        }
        if (a.empty() || b.empty()) {
            return a.size() > b.size();
        }
        return a[0] < b[0];
    });
}

static void loo_scores(const std::vector<int>& members, const std::vector<Pose>& poses, Method method,
                       std::vector<float>* scores) {
    int n = static_cast<int>(members.size());
    scores->assign(static_cast<size_t>(n), 0.f);
    int threads = thread_count(n);
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t) {
        pool.emplace_back([&, t]() {
            Work work;
            prepare_work(&work);
            std::vector<int> protos;
            protos.reserve(static_cast<size_t>(n));
            for (int i = t; i < n; i += threads) {
                protos.clear();
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        protos.push_back(members[static_cast<size_t>(j)]);
                    }
                }
                (*scores)[static_cast<size_t>(i)] = score_query(members[static_cast<size_t>(i)], protos, poses, method, &work);
            }
        });
    }
    for (size_t t = 0; t < pool.size(); ++t) {
        pool[t].join();
    }
}

static float median_score(std::vector<float> scores) {
    asserta(!scores.empty());
    std::sort(scores.begin(), scores.end());
    int n = static_cast<int>(scores.size());
    if ((n % 2) == 1) {
        return scores[static_cast<size_t>(n / 2)];
    }
    return 0.5f * (scores[static_cast<size_t>(n / 2 - 1)] + scores[static_cast<size_t>(n / 2)]);
}

static void report_small_clusters(const std::vector<std::vector<int>>& comps, const std::vector<char>& is_main,
                                  int decoy_begin, int n_decoy, Method method, int trans, int rot180, int scale,
                                  const char* path, int* scored, int* singletons, int* union_hits) {
    std::vector<Pose> poses;
    make_poses(trans, rot180, scale, &poses);
    std::vector<int> decoys(static_cast<size_t>(n_decoy));
    for (int i = 0; i < n_decoy; ++i) {
        decoys[static_cast<size_t>(i)] = decoy_begin + i;
    }
    std::vector<char> hit(static_cast<size_t>(n_decoy), 0);
    FILE* out = std::fopen(path, "wb");
    if (!out) {
        die("cannot write small_clusters.tsv");
    }
    std::fprintf(out,
                 "rank\tsize\tscored\tloo_max\tloo_median\tdecoys_covered\tdecoy_n\tdecoy_rate\t"
                 "members_at_1pct\n");
    std::fflush(out);
    *scored = 0;
    *singletons = 0;
    int rank = 0;
    for (size_t c = 0; c < comps.size(); ++c) {
        const std::vector<int>& members = comps[c];
        if (members.empty() || is_main[static_cast<size_t>(members[0])] != 0) {
            continue;
        }
        rank += 1;
        int size = static_cast<int>(members.size());
        if (size < 2) {
            *singletons += 1;
            std::fprintf(out, "%d\t%d\t0\t-\t-\t-\t%d\t-\t-\n", rank, size, n_decoy);
            std::fprintf(stderr, "small cluster rank %d size %d singleton\n", rank, size);
            std::fflush(out);
            continue;
        }
        std::fprintf(stderr, "small cluster rank %d size %d\n", rank, size);
        std::vector<float> member_s;
        loo_scores(members, poses, method, &member_s);
        float loo_max = member_s[0];
        for (size_t i = 1; i < member_s.size(); ++i) {
            if (member_s[i] > loo_max) {
                loo_max = member_s[i];
            }
        }
        float loo_med = median_score(member_s);
        std::vector<float> decoy_s;
        score_list(decoys, members, poses, method, &decoy_s);
        int covered = count_le(decoy_s, loo_max);
        float fpr1 = threshold_at_fpr(decoy_s, 0.01);
        int members_at_1 = count_le(member_s, fpr1);
        double rate = static_cast<double>(covered) / static_cast<double>(n_decoy);
        for (int i = 0; i < n_decoy; ++i) {
            if (decoy_s[static_cast<size_t>(i)] <= loo_max) {
                hit[static_cast<size_t>(i)] = 1;
            }
        }
        *scored += 1;
        std::fprintf(out, "%d\t%d\t1\t%.6g\t%.6g\t%d\t%d\t%.6f\t%d\n", rank, size, loo_max, loo_med, covered,
                     n_decoy, rate, members_at_1);
        std::fprintf(stderr, "small cluster rank %d loo_max %.4f decoys %d/%d members_at_1pct %d/%d\n", rank,
                     loo_max, covered, n_decoy, members_at_1, size);
        std::fflush(out);
    }
    std::fclose(out);
    int unions = 0;
    for (int i = 0; i < n_decoy; ++i) {
        if (hit[static_cast<size_t>(i)] != 0) {
            unions += 1;
        }
    }
    *union_hits = unions;
}

static void ensure_dir(const char* path) {
#if defined(_WIN32)
    _mkdir(path);
#else
    mkdir(path, 0755);
#endif
}

static void usage() {
    std::fprintf(stdout,
                 "usage: avsun_classifier --avsun DIR --decoy DIR --out DIR\n"
                 "\n"
                 "Search one-class thumbnail designs on the largest avsun shape cluster.\n"
                 "A low score matches a prototype from that cluster. The chosen design\n"
                 "has the highest validation recall at a 1%% decoy false-positive rate.\n");
}

int main(int argc, char** argv) {
    const char* avsun_dir = 0;
    const char* decoy_dir = 0;
    const char* out_dir = 0;
    for (int i = 1; i < argc; ++i) {
        const char* arg = argv[i];
        if (std::strcmp(arg, "--help") == 0) {
            usage();
            return 0;
        }
        const char** slot = 0;
        if (std::strcmp(arg, "--avsun") == 0) {
            slot = &avsun_dir;
        } else if (std::strcmp(arg, "--decoy") == 0) {
            slot = &decoy_dir;
        } else if (std::strcmp(arg, "--out") == 0) {
            slot = &out_dir;
        } else {
            std::fprintf(stderr, "avsun_classifier: unknown argument %s\n", arg);
            return 1;
        }
        if (i + 1 >= argc) {
            die("missing option value");
        }
        *slot = argv[++i];
    }
    if (!avsun_dir || !decoy_dir || !out_dir) {
        usage();
        return 1;
    }
    self_test_edt();
    int avsun_begin = 0;
    load_dir(avsun_dir, &avsun_begin);
    int n_avsun = static_cast<int>(g_im.size());
    int decoy_begin = 0;
    load_dir(decoy_dir, &decoy_begin);
    int n_decoy = static_cast<int>(g_im.size()) - decoy_begin;
    if (n_avsun < 4 || n_decoy < 4) {
        die("need at least 4 avsun and 4 decoy thumbnails");
    }
    std::vector<int> all_ids(static_cast<size_t>(n_avsun));
    for (int i = 0; i < n_avsun; ++i) {
        all_ids[static_cast<size_t>(i)] = i;
    }
    std::vector<Link> links;
    collect_links(all_ids, kLink, &links);
    std::vector<int> group_of;
    int n_groups = 0;
    int largest = 0;
    groups_at(n_avsun, links, kLink, &group_of, &n_groups, &largest);
    std::fprintf(stderr, "clusters %d\n", n_groups);
    std::fprintf(stderr, "largest cluster %d of %d\n", largest, n_avsun);
    std::vector<std::vector<int>> shape_clusters;
    list_components(group_of, &shape_clusters);
    std::vector<int> kept;
    std::vector<int> dropped;
    largest_component(group_of, &kept, &dropped);
    std::vector<char> is_main(static_cast<size_t>(n_avsun), 0);
    for (size_t i = 0; i < kept.size(); ++i) {
        is_main[static_cast<size_t>(kept[i])] = 1;
    }
    links.clear();
    const float kCuts[] = {0.f, 0.05f, 0.10f, 0.25f, 0.50f, 0.75f, 1.00f, 1.25f, 1.50f, 1.75f};
    const int n_cuts = static_cast<int>(sizeof(kCuts) / sizeof(kCuts[0]));
    collect_links(kept, kCuts[n_cuts - 1], &links);
    int kept_n = static_cast<int>(kept.size());
    int train_cap = (kept_n + 1) / 2;
    std::vector<int> train;
    std::vector<int> val_av;
    std::vector<int> test_av;
    float chosen_cut = -1.f;
    int chosen_groups = 0;
    int chosen_largest = 0;
    for (int c = 0; c < n_cuts; ++c) {
        int sub_groups = 0;
        int sub_largest = 0;
        groups_at(kept_n, links, kCuts[c], &group_of, &sub_groups, &sub_largest);
        std::fprintf(stderr, "subcluster cut %.2f groups %d largest %d\n", kCuts[c], sub_groups, sub_largest);
        if (sub_largest > train_cap) {
            break;
        }
        std::vector<int> try_train;
        std::vector<int> try_val;
        std::vector<int> try_test;
        if (!split_groups(kept, group_of, &try_train, &try_val, &try_test)) {
            std::fprintf(stderr, "subcluster cut %.2f leaves an empty fold\n", kCuts[c]);
            break;
        }
        chosen_cut = kCuts[c];
        chosen_groups = sub_groups;
        chosen_largest = sub_largest;
        train.swap(try_train);
        val_av.swap(try_val);
        test_av.swap(try_test);
    }
    if (chosen_cut < 0.f) {
        die("tight subclusters still leave an empty avsun fold");
    }
    std::fprintf(stderr, "subcluster chosen cut %.2f groups %d largest %d\n", chosen_cut, chosen_groups,
                 chosen_largest);
    std::fprintf(stderr, "split avsun train %zu val %zu test %zu discarded %zu\n", train.size(), val_av.size(),
                 test_av.size(), dropped.size());
    float radius = medoid_radius(kept);
    std::fprintf(stderr, "medoid radius %.4f\n", radius);
    std::vector<int> val_de;
    std::vector<int> test_de;
    split_decoys(decoy_begin, n_decoy, &val_de, &test_de);

    ensure_dir(out_dir);
    std::string results_path = std::string(out_dir) + "/results.tsv";
    std::string summary_path = std::string(out_dir) + "/summary.txt";
    FILE* results = std::fopen(results_path.c_str(), "wb");
    if (!results) {
        die("cannot write results.tsv");
    }
    std::fprintf(results,
                 "method\tk\ttrans\trot180\tscale\tval_n_avsun\tval_n_decoy\tval_recall_fpr1\tval_fp\t"
                 "test_n_avsun\ttest_n_decoy\ttest_recall_fpr1\ttest_fp\ttest_recall_fpr5\t"
                 "test_recall_fpr01\tthreshold\n");
    std::fflush(results);

    const int ks[3] = {1, 8, 32};
    Design best;
    best.val_recall_fpr1 = -1.0;
    best.method = kChamfer;
    best.k = 1;
    best.trans = 0;
    best.rot180 = 0;
    best.scale = 0;
    best.val_fp = 0;
    best.test_recall_fpr1 = 0;
    best.test_fp = 0;
    best.test_recall_fpr5 = 0;
    best.test_recall_fpr01 = 0;
    best.threshold = 0;
    best.val_hits = 0;
    best.test_hits = 0;
    bool have_best = false;

    for (int mi = 0; mi < 3; ++mi) {
        Method method = static_cast<Method>(mi);
        for (int ki = 0; ki < 3; ++ki) {
            int k = ks[ki];
            std::fprintf(stderr, "prototypes %s k=%d\n", method_name(method), k);
            std::vector<int> protos = farthest_prototypes(method, train, k);
            for (int trans = 0; trans <= 1; ++trans) {
                for (int rot180 = 0; rot180 <= 1; ++rot180) {
                    for (int scale = 0; scale <= 1; ++scale) {
                        std::vector<Pose> poses;
                        make_poses(trans, rot180, scale, &poses);
                        std::fprintf(stderr, "score %s k=%d trans=%d rot180=%d scale=%d poses=%zu\n",
                                     method_name(method), static_cast<int>(protos.size()), trans, rot180,
                                     scale, poses.size());
                        std::vector<float> val_av_s;
                        std::vector<float> val_de_s;
                        std::vector<float> test_av_s;
                        std::vector<float> test_de_s;
                        score_list(val_av, protos, poses, method, &val_av_s);
                        score_list(val_de, protos, poses, method, &val_de_s);
                        score_list(test_av, protos, poses, method, &test_av_s);
                        score_list(test_de, protos, poses, method, &test_de_s);
                        float thr1 = threshold_at_fpr(val_de_s, 0.01);
                        float test_thr5 = threshold_at_fpr(test_de_s, 0.05);
                        float test_thr01 = threshold_at_fpr(test_de_s, 0.001);
                        Design row;
                        row.method = method;
                        row.k = static_cast<int>(protos.size());
                        row.trans = trans;
                        row.rot180 = rot180;
                        row.scale = scale;
                        row.threshold = thr1;
                        row.val_recall_fpr1 = recall_le(val_av_s, thr1);
                        row.val_hits = count_le(val_av_s, thr1);
                        row.val_fp = count_le(val_de_s, thr1);
                        row.test_recall_fpr1 = recall_le(test_av_s, thr1);
                        row.test_hits = count_le(test_av_s, thr1);
                        row.test_fp = count_le(test_de_s, thr1);
                        row.test_recall_fpr5 = recall_le(test_av_s, test_thr5);
                        row.test_recall_fpr01 = recall_le(test_av_s, test_thr01);
                        std::fprintf(results,
                                     "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%d\t%d\t%d\t%.6f\t%d\t%.6f\t%.6f\t"
                                     "%.6g\n",
                                     method_name(method), row.k, trans, rot180, scale,
                                     static_cast<int>(val_av.size()), static_cast<int>(val_de.size()),
                                     row.val_recall_fpr1, row.val_fp, static_cast<int>(test_av.size()),
                                     static_cast<int>(test_de.size()), row.test_recall_fpr1, row.test_fp,
                                     row.test_recall_fpr5, row.test_recall_fpr01, row.threshold);
                        std::fflush(results);
                        bool better = !have_best || row.val_recall_fpr1 > best.val_recall_fpr1 + 1.0e-12;
                        if (!better && have_best &&
                            std::fabs(row.val_recall_fpr1 - best.val_recall_fpr1) <= 1.0e-12) {
                            int row_pose = trans + rot180 + scale;
                            int best_pose = best.trans + best.rot180 + best.scale;
                            if (row.k < best.k || (row.k == best.k && row_pose < best_pose) ||
                                (row.k == best.k && row_pose == best_pose && row.method < best.method)) {
                                better = true;
                            }
                        }
                        if (better) {
                            best = row;
                            have_best = true;
                        }
                    }
                }
            }
        }
    }
    std::fclose(results);
    std::vector<int> win_protos = farthest_prototypes(best.method, train, best.k);
    std::vector<Pose> win_poses;
    make_poses(best.trans, best.rot180, best.scale, &win_poses);
    std::vector<float> dropped_s;
    score_list(dropped, win_protos, win_poses, best.method, &dropped_s);
    int dropped_hits = count_le(dropped_s, best.threshold);
    std::fprintf(stderr, "discarded passing %d of %zu\n", dropped_hits, dropped.size());
    int small_scored = 0;
    int small_singletons = 0;
    int small_union = 0;
    std::string small_path = std::string(out_dir) + "/small_clusters.tsv";
    report_small_clusters(shape_clusters, is_main, decoy_begin, n_decoy, best.method, best.trans, best.rot180,
                          best.scale, small_path.c_str(), &small_scored, &small_singletons, &small_union);
    FILE* summary = std::fopen(summary_path.c_str(), "wb");
    if (!summary) {
        die("cannot write summary.txt");
    }
    std::fprintf(summary, "winner %s k=%d trans=%d rot180=%d scale=%d\n", method_name(best.method), best.k,
                 best.trans, best.rot180, best.scale);
    std::fprintf(summary, "positive avsuns %d of %d at chamfer <= %.0f\n", kept_n, n_avsun, kLink);
    std::fprintf(summary, "subclusters at chamfer <= %.2f: %d groups, largest %d\n", chosen_cut, chosen_groups,
                 chosen_largest);
    std::fprintf(summary, "avsun train %d val %d test %d\n", static_cast<int>(train.size()),
                 static_cast<int>(val_av.size()), static_cast<int>(test_av.size()));
    std::fprintf(summary, "discarded avsuns %d\n", static_cast<int>(dropped.size()));
    std::fprintf(summary, "medoid radius %.4f\n", radius);
    std::fprintf(summary, "decoy val %d test %d\n", static_cast<int>(val_de.size()),
                 static_cast<int>(test_de.size()));
    std::fprintf(summary, "threshold %.6g\n", best.threshold);
    std::fprintf(summary, "val recall at 1%% FPR %.6f (%d of %d avsuns, %d false positives of %d decoys)\n",
                 best.val_recall_fpr1, best.val_hits, static_cast<int>(val_av.size()), best.val_fp,
                 static_cast<int>(val_de.size()));
    std::fprintf(summary,
                 "test recall at that threshold %.6f (%d of %d avsuns, %d false positives of %d decoys)\n",
                 best.test_recall_fpr1, best.test_hits, static_cast<int>(test_av.size()), best.test_fp,
                 static_cast<int>(test_de.size()));
    std::fprintf(summary, "test recall at 5%% FPR %.6f\n", best.test_recall_fpr5);
    std::fprintf(summary, "test recall at 0.1%% FPR %.6f\n", best.test_recall_fpr01);
    std::fprintf(summary, "discarded avsuns passing threshold %d of %d\n", dropped_hits,
                 static_cast<int>(dropped.size()));
    std::fprintf(summary, "small clusters scored %d singletons %d\n", small_scored, small_singletons);
    std::fprintf(summary, "small-cluster union %d of %d decoys (%.6f)\n", small_union, n_decoy,
                 n_decoy > 0 ? static_cast<double>(small_union) / static_cast<double>(n_decoy) : 0.0);
    std::fclose(summary);
    std::printf("winner %s k=%d trans=%d rot180=%d scale=%d\n", method_name(best.method), best.k, best.trans,
                best.rot180, best.scale);
    std::printf("val recall at 1%% FPR %.6f\n", best.val_recall_fpr1);
    std::printf("test recall %.6f (%d of %d avsuns) with %d / %d decoys passing\n", best.test_recall_fpr1,
                best.test_hits, static_cast<int>(test_av.size()), best.test_fp, static_cast<int>(test_de.size()));
    std::printf("discarded passing %d of %d\n", dropped_hits, static_cast<int>(dropped.size()));
    std::printf("small-cluster union %d of %d decoys\n", small_union, n_decoy);
    std::printf("wrote %s\n", results_path.c_str());
    std::printf("wrote %s\n", summary_path.c_str());
    std::printf("wrote %s\n", small_path.c_str());
    return 0;
}
