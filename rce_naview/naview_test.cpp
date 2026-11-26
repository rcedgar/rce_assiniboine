#include "filename.hpp"
#include "naview.hpp"
#include "svg_page.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

static int g_fails = 0;

static void fail(const char* msg) {
    std::fprintf(stderr, "FAIL %s\n", msg);
    g_fails += 1;
}

static void expect(bool cond, const char* msg) {
    if (!cond) {
        fail(msg);
    }
}

static void expect_str(const std::string& got, const char* want, const char* msg) {
    if (got != want) {
        std::fprintf(stderr, "FAIL %s\n  got  %s\n  want %s\n", msg, got.c_str(), want);
        g_fails += 1;
    }
}

static bool near(double a, double b) {
    return std::fabs(a - b) <= 1e-6;
}

static bool has_line(const NavPicture& pic, double x0, double y0, double x1, double y1) {
    for (size_t i = 0; i < pic.lines.size(); ++i) {
        const NavLine& L = pic.lines[i];
        if (near(L.x0, x0) && near(L.y0, y0) && near(L.x1, x1) && near(L.y1, y1)) {
            return true;
        }
        if (near(L.x0, x1) && near(L.y0, y1) && near(L.x1, x0) && near(L.y1, y0)) {
            return true;
        }
    }
    return false;
}

static void expect_line(const NavPicture& pic, double x0, double y0, double x1, double y1, const char* msg) {
    if (!has_line(pic, x0, y0, x1, y1)) {
        std::fprintf(stderr, "FAIL %s missing line (%.6g,%.6g)-(%.6g,%.6g) among %zu\n", msg, x0, y0,
                     x1, y1, pic.lines.size());
        for (size_t i = 0; i < pic.lines.size() && i < 12; ++i) {
            const NavLine& L = pic.lines[i];
            std::fprintf(stderr, "  (%.6g,%.6g)-(%.6g,%.6g)\n", L.x0, L.y0, L.x1, L.y1);
        }
        g_fails += 1;
    }
}

static NavLayout must_layout(const char* db, bool circular) {
    NavLayout lay;
    const char* err = 0;
    bool ok = naview_layout(db, static_cast<int>(std::strlen(db)), circular, &lay, &err);
    if (!ok) {
        std::fprintf(stderr, "FAIL layout %s circular=%d err=%s\n", db, circular ? 1 : 0,
                     err ? err : "");
        g_fails += 1;
    }
    return lay;
}

static void expect_reject(const char* db) {
    NavLayout lay;
    const char* err = 0;
    bool ok = naview_layout(db, static_cast<int>(std::strlen(db)), false, &lay, &err);
    if (ok) {
        std::fprintf(stderr, "FAIL expected reject for %s\n", db);
        g_fails += 1;
    }
}

static void check_pairs(const NavLayout& lay, const char* msg) {
    for (size_t i = 0; i < lay.mate.size(); ++i) {
        int j = lay.mate[i];
        if (j < 0 || j < static_cast<int>(i)) {
            continue;
        }
        double dx = lay.x[i] - lay.x[static_cast<size_t>(j)];
        double dy = lay.y[i] - lay.y[static_cast<size_t>(j)];
        double d = std::sqrt(dx * dx + dy * dy);
        if (std::fabs(d - 1.0) > 1e-6) {
            std::fprintf(stderr, "FAIL %s pair %zu-%d distance %.8g\n", msg, i, j, d);
            g_fails += 1;
            return;
        }
    }
}

static void test_worked_linear() {
    const char* db = ".((..((...))..((...)).))";
    NavLayout lay = must_layout(db, false);
    if (lay.loops.size() != 4 || lay.helices.size() != 3) {
        fail("worked linear counts");
        return;
    }
    expect(lay.loops[0].kind == NavKind::Exterior, "L0 kind");
    expect(lay.loops[0].start == 0 && lay.loops[0].end == 23, "L0 span");
    expect(lay.loops[0].parent < 0 && lay.loops[0].closing_helix < 0, "L0 links");
    expect(lay.loops[0].child_helices.size() == 1 && lay.loops[0].child_helices[0] == 0, "L0 child");

    expect(lay.loops[1].kind == NavKind::Multibranch, "L1 kind");
    expect(lay.loops[1].start == 2 && lay.loops[1].end == 22, "L1 span");
    expect(lay.loops[1].parent == 0 && lay.loops[1].closing_helix == 0, "L1 links");
    expect(lay.loops[1].child_helices.size() == 2, "L1 children");
    expect(lay.loops[1].child_helices[0] == 1 && lay.loops[1].child_helices[1] == 2, "L1 child ids");

    expect(lay.loops[2].kind == NavKind::Hairpin && lay.loops[2].start == 6 && lay.loops[2].end == 10,
           "L2");
    expect(lay.loops[3].kind == NavKind::Hairpin && lay.loops[3].start == 15 && lay.loops[3].end == 19,
           "L3");

    expect(lay.helices[0].start5 == 1 && lay.helices[0].end5 == 2, "H0 5p");
    expect(lay.helices[0].start3 == 22 && lay.helices[0].end3 == 23, "H0 3p");
    expect(lay.helices[1].start5 == 5 && lay.helices[1].end5 == 6, "H1 5p");
    expect(lay.helices[1].start3 == 10 && lay.helices[1].end3 == 11, "H1 3p");
    expect(lay.helices[2].start5 == 14 && lay.helices[2].end5 == 15, "H2 5p");
    expect(lay.helices[2].start3 == 19 && lay.helices[2].end3 == 20, "H2 3p");

    struct Want {
        NavKind kind;
        int start;
        int end;
        int helix;
        int loop;
        NavSide side;
    };
    const Want want[] = {
        {NavKind::Exterior, 0, 0, -1, 0, NavSide::None},
        {NavKind::Helix5, 1, 2, 0, -1, NavSide::Five},
        {NavKind::Multibranch, 3, 4, -1, 1, NavSide::None},
        {NavKind::Helix5, 5, 6, 1, -1, NavSide::Five},
        {NavKind::Hairpin, 7, 9, -1, 2, NavSide::None},
        {NavKind::Helix3, 10, 11, 1, -1, NavSide::Three},
        {NavKind::Multibranch, 12, 13, -1, 1, NavSide::None},
        {NavKind::Helix5, 14, 15, 2, -1, NavSide::Five},
        {NavKind::Hairpin, 16, 18, -1, 3, NavSide::None},
        {NavKind::Helix3, 19, 20, 2, -1, NavSide::Three},
        {NavKind::Multibranch, 21, 21, -1, 1, NavSide::None},
        {NavKind::Helix3, 22, 23, 0, -1, NavSide::Three},
    };
    expect(lay.segments.size() == 12, "segment count");
    if (lay.segments.size() == 12) {
        for (int i = 0; i < 12; ++i) {
            const NavSegment& s = lay.segments[static_cast<size_t>(i)];
            if (s.kind != want[i].kind || s.start != want[i].start || s.end != want[i].end ||
                s.helix_id != want[i].helix || s.loop_id != want[i].loop || s.side != want[i].side) {
                std::fprintf(stderr,
                             "FAIL segment %d got kind %d %d-%d h %d l %d side %d\n", i,
                             static_cast<int>(s.kind), s.start, s.end, s.helix_id, s.loop_id,
                             static_cast<int>(s.side));
                g_fails += 1;
            }
        }
    }
    expect_str(lay.naview_graph, "E1(2:M2+2+1(2:P3,2:P3))", "naview graph linear");
    expect_str(lay.rag_tree, "(*,@,@)", "rag tree linear");
    expect_str(lay.rag_dual, "(*-@-@-)", "rag dual linear");

    const double s3 = std::sqrt(3.0);
    expect_line(lay.rag_picture, 0, 0, 0.5, s3 / 2.0, "rag tree star");
    expect_line(lay.rag_picture, 0, 0, -1, 0, "rag tree first hairpin");
    expect_line(lay.rag_picture, 0, 0, 0.5, -s3 / 2.0, "rag tree second hairpin");
    expect_line(lay.naview_picture, 0, 0, -1, 0, "naview graph stem");
    expect_line(lay.naview_picture, -1, 0, -1, 1, "naview graph hairpin a");
    expect_line(lay.naview_picture, -1, 0, -1, -1, "naview graph hairpin b");
    expect_line(lay.dual_picture, 1, 0, -0.5, s3 / 2.0, "dual 0-1");
    expect_line(lay.dual_picture, -0.5, s3 / 2.0, -0.5, -s3 / 2.0, "dual 1-2");
    expect_line(lay.dual_picture, -0.5, -s3 / 2.0, 1, 0, "dual 2-0");
    expect_line(lay.dual_picture, -0.725, 0.725 * s3, (-5.0 - s3) / 8.0, (5.0 * s3 - 1.0) / 8.0,
                "dual triangle");
    check_pairs(lay, "worked linear");
}

static void test_worked_circular() {
    const char* db = ".((..((...))..((...)).))";
    NavLayout lay = must_layout(db, true);
    expect(lay.loops[0].kind == NavKind::Hairpin, "circ kind");
    expect(lay.loops[0].start == 0 && lay.loops[0].end == 0, "circ span");
    expect(lay.loops[0].parent < 0 && lay.loops[0].closing_helix < 0, "circ links");
    expect_str(lay.naview_graph, "P1(2:M2+2+1(2:P3,2:P3))", "naview graph circ");
    expect_str(lay.rag_tree, "(@,@)", "rag tree circ");
    expect_str(lay.rag_dual, "(@=@)", "rag dual circ");
    expect(lay.segments[0].kind == NavKind::Hairpin && lay.segments[0].start == 0 &&
               lay.segments[0].end == 0,
           "circ base0 hairpin");
    bool exterior = false;
    for (size_t i = 0; i < lay.segments.size(); ++i) {
        if (lay.segments[i].kind == NavKind::Exterior) {
            exterior = true;
        }
    }
    expect(!exterior, "circ has no exterior segment");
    expect_line(lay.rag_picture, 0, 0, 0, 1, "circ rag a");
    expect_line(lay.rag_picture, 0, 0, 0, -1, "circ rag b");
    expect_line(lay.dual_picture, 1, -0.2, -1, -0.2, "circ dual +p");
    expect_line(lay.dual_picture, 1, 0.2, -1, 0.2, "circ dual -p");
    expect_line(lay.dual_picture, 1.45, 0, 1.25, 0.25, "circ triangle");
    expect_line(lay.dual_picture, -1.45, 0, -1.25, -0.25, "circ triangle other");
    check_pairs(lay, "worked circular");
}

static void test_interior_bulge() {
    NavLayout interior = must_layout("((..((...))..))", false);
    expect_str(interior.naview_graph, "E0(2:I2+2(2:P3))", "interior graph");
    expect_str(interior.rag_tree, "(*,@)", "interior rag");
    expect_str(interior.rag_dual, "(*=@)", "interior dual");
    bool saw5 = false;
    bool saw3 = false;
    for (size_t i = 0; i < interior.segments.size(); ++i) {
        const NavSegment& s = interior.segments[i];
        if (s.kind == NavKind::Interior && s.side == NavSide::Five && s.start == 2 && s.end == 3) {
            saw5 = true;
        }
        if (s.kind == NavKind::Interior && s.side == NavSide::Three && s.start == 11 && s.end == 12) {
            saw3 = true;
        }
    }
    expect(saw5 && saw3, "interior sides");

    NavLayout bulge = must_layout("((..((...))))", false);
    expect_str(bulge.naview_graph, "E0(2:B2+0(2:P3))", "bulge graph");
    expect_str(bulge.rag_tree, "(*,@)", "bulge rag");
    expect_str(bulge.rag_dual, "(*-@)", "bulge dual");
    int bulge_segs = 0;
    bool bulge_3 = false;
    for (size_t i = 0; i < bulge.segments.size(); ++i) {
        const NavSegment& s = bulge.segments[i];
        if (s.kind == NavKind::Bulge) {
            bulge_segs += 1;
            if (s.side == NavSide::Five && s.start == 2 && s.end == 3) {
                // expected sole bulge segment
            } else {
                fail("unexpected bulge segment");
            }
            if (s.side == NavSide::Three) {
                bulge_3 = true;
            }
        }
    }
    expect(bulge_segs == 1 && !bulge_3, "bulge omits empty side");
    check_pairs(interior, "interior");
    check_pairs(bulge, "bulge");
}

static void test_no_pairs() {
    NavLayout one = must_layout(".", false);
    expect(one.loops.size() == 1 && one.loops[0].kind == NavKind::Exterior, "n1 kind");
    expect(one.x.size() == 1 && near(one.x[0], 0) && near(one.y[0], 0), "n1 origin");
    expect_str(one.naview_graph, "E1", "n1 graph");
    expect_str(one.rag_tree, "*", "n1 rag");
    expect_str(one.rag_dual, "*", "n1 dual");
    expect(one.dual_picture.vertices.size() == 1 && one.dual_picture.lines.empty(), "lone star");

    NavLayout circ = must_layout("....", true);
    expect(circ.loops[0].kind == NavKind::Circle, "circle kind");
    expect_str(circ.naview_graph, "C4", "circle graph");
    expect_str(circ.rag_tree, "*", "circle rag");
    expect_str(circ.rag_dual, "*", "circle dual");
    expect(circ.segments.size() == 1 && circ.segments[0].kind == NavKind::Circle, "circle seg");
    for (int i = 0; i < 4; ++i) {
        int j = (i + 1) % 4;
        double dx = circ.x[static_cast<size_t>(i)] - circ.x[static_cast<size_t>(j)];
        double dy = circ.y[static_cast<size_t>(i)] - circ.y[static_cast<size_t>(j)];
        double d = std::sqrt(dx * dx + dy * dy);
        if (std::fabs(d - 1.0) > 1e-6) {
            std::fprintf(stderr, "FAIL circle chord %.8g\n", d);
            g_fails += 1;
            break;
        }
    }
    NavLayout lin = must_layout("....", false);
    expect_str(lin.naview_graph, "E4", "linear unpaired graph");
}

static void test_reject() {
    expect_reject("(");
    expect_reject("())");
    expect_reject("(]");
    expect_reject("((())");
}

static std::string read_all(const char* path) {
    std::string data;
    FILE* f = std::fopen(path, "rb");
    if (!f) {
        return data;
    }
    char buf[4096];
    size_t n = 0;
    while ((n = std::fread(buf, 1, sizeof(buf), f)) > 0) {
        data.append(buf, n);
    }
    std::fclose(f);
    return data;
}

static bool file_contains(const char* path, const char* needle) {
    return read_all(path).find(needle) != std::string::npos;
}

static bool backbone_span(const std::string& data, double* dx, double* dy, double* width) {
    size_t at = 0;
    while ((at = data.find("<line ", at)) != std::string::npos) {
        size_t end = data.find("/>", at);
        if (end == std::string::npos) {
            return false;
        }
        std::string line = data.substr(at, end - at);
        if (line.find("stroke=\"#212121\"") != std::string::npos) {
            double x1 = 0.0;
            double y1 = 0.0;
            double x2 = 0.0;
            double y2 = 0.0;
            if (std::sscanf(line.c_str(),
                            "<line x1=\"%lf\" y1=\"%lf\" x2=\"%lf\" y2=\"%lf\" "
                            "stroke=\"#212121\" stroke-width=\"%lf\"",
                            &x1, &y1, &x2, &y2, width) != 5) {
                return false;
            }
            *dx = x2 - x1;
            *dy = y2 - y1;
            return true;
        }
        at += 1;
    }
    return false;
}

static void test_svg_smoke() {
#if defined(_WIN32)
    _mkdir("build-linux");
    _mkdir("build-linux/naview_test_out");
#else
    mkdir("build-linux", 0755);
    mkdir("build-linux/naview_test_out", 0755);
#endif
    const char* db = "(((...)))";
    const char* seq = "GGGAAACCC";
    NavLayout lay = must_layout(db, false);
    const char* path = "build-linux/naview_test_out/smoke.svg";
    expect(write_naview_svg(path, "A&B", seq, db, false, &lay, 60, 1.0, false, false, 0, 0, 0),
           "svg write");
    const char* only_path = "build-linux/naview_test_out/only.svg";
    expect(write_naview_svg(only_path, "A&B", seq, db, false, &lay, 60, 1.0, true, false, 0, 0, 0),
           "naview only");
    expect(!file_contains(only_path, "id=\"title\""), "only omits title");
    expect(!file_contains(only_path, "id=\"dotbracket\""), "only omits dot-bracket");
    expect(!file_contains(only_path, "id=\"graph-naview\""), "only omits graphs");
    expect(file_contains(only_path, "id=\"naview\""), "only keeps naview");
    const char* fit_path = "build-linux/naview_test_out/fit.svg";
    const char* fit_only_path = "build-linux/naview_test_out/fit_only.svg";
    expect(write_naview_svg(fit_path, "A&B", seq, db, false, &lay, 2, 1.0, false, false, 0, 0, 0),
           "fitted write");
    expect(write_naview_svg(fit_only_path, "A&B", seq, db, false, &lay, 2, 1.0, true, false, 0, 0, 0),
           "fitted only");
    double fit_dx = 0.0;
    double fit_dy = 0.0;
    double fit_w = 0.0;
    double only_dx = 0.0;
    double only_dy = 0.0;
    double only_w = 0.0;
    expect(backbone_span(read_all(fit_path), &fit_dx, &fit_dy, &fit_w), "fitted backbone");
    expect(backbone_span(read_all(fit_only_path), &only_dx, &only_dy, &only_w), "only backbone");
    expect(near(fit_dx, only_dx) && near(fit_dy, only_dy) && near(fit_w, only_w),
           "only matches fitted naview");
    expect(fit_w < 1.25, "naview scaled below natural stroke");
    expect(file_contains(path, "id=\"title\""), "title group");
    expect(file_contains(path, "A&amp;B"), "escaped title");
    expect(file_contains(path, "id=\"naview\""), "naview group");
    expect(file_contains(path, "id=\"dotbracket\""), "db group");
    expect(file_contains(path, "id=\"graph-naview\""), "graph naview");
    expect(file_contains(path, "id=\"graph-rag\""), "graph rag");
    expect(file_contains(path, "id=\"graph-rag-dual\""), "graph dual");
    expect(file_contains(path, "textLength="), "forced text width");
    expect(file_contains(path, "fill-opacity=\"0.10\""), "wash opacity");
    expect(file_contains(path, "ui-monospace"), "monospace");
    expect(file_contains(path, "#c62828"), "default red");
    expect(file_contains(path, "#1565c0"), "default blue");
    const char* rain = "build-linux/naview_test_out/rainbow.svg";
    expect(write_naview_svg(rain, "A&B", seq, db, false, &lay, 60, 1.0, false, true, 0, 0, 0),
           "rainbow write");
    expect(file_contains(rain, "hsl(0.00, 72%, 40%)"), "rainbow first hue");
    const char* custom_cols[] = {"#112233", "lime"};
    const char* custom = "build-linux/naview_test_out/custom.svg";
    expect(write_naview_svg(custom, "A&B", seq, db, false, &lay, 60, 1.0, false, false, custom_cols,
                            2, 0),
           "custom write");
    expect(file_contains(custom, "#112233"), "custom first");
    expect(file_contains(custom, "lime"), "custom second");
    check_pairs(lay, "smoke hairpin");
    const char* window_path = "build-linux/naview_test_out/window.svg";
    NavLayout guanine = must_layout(".........", false);
    expect(write_naview_svg(window_path, "G", "GGGGGGGGG", ".........", false, &guanine, 60, 1.0,
                            true, false, 0, 0, 20),
           "window write");
    expect(file_contains(window_path, "#ffa360"), "all-G composition");
    expect(!file_contains(window_path, "hsl("), "window is not rainbow");
}

static void test_svg_color() {
    expect(naview_svg_color("red"), "name");
    expect(naview_svg_color("RED"), "name case");
    expect(naview_svg_color("  blue  "), "trimmed name");
    expect(naview_svg_color("transparent"), "transparent");
    expect(naview_svg_color("currentColor"), "currentcolor");
    expect(naview_svg_color("#abc"), "hex3");
    expect(naview_svg_color("#abcd"), "hex4");
    expect(naview_svg_color("#aabbcc"), "hex6");
    expect(naview_svg_color("#aabbccdd"), "hex8");
    expect(naview_svg_color("rgb(1, 2, 3)"), "rgb");
    expect(naview_svg_color("rgb(1e2, 0, 0)"), "rgb exponent");
    expect(naview_svg_color("rgb(10%, 0%, 100%)"), "rgb percent");
    expect(naview_svg_color("rgba(1,2,3,0.5)"), "rgba");
    expect(naview_svg_color("hsl(120, 50%, 40%)"), "hsl");
    expect(naview_svg_color("hsla(120, 50%, 40%, 50%)"), "hsla");
    expect(!naview_svg_color(""), "empty");
    expect(!naview_svg_color("notacolor"), "unknown name");
    expect(!naview_svg_color("#gg0000"), "bad hex");
    expect(!naview_svg_color("#ab"), "short hex");
    expect(!naview_svg_color("rgb(1,2)"), "rgb arity");
    expect(!naview_svg_color("rgb(1, 2%, 3)"), "mixed percent");
    expect(!naview_svg_color("hsl(120%, 50%, 40%)"), "hsl hue percent");
}

static void test_wrap_segment() {
    NavLayout lay = must_layout(".(...).", true);
    bool found = false;
    for (size_t i = 0; i < lay.segments.size(); ++i) {
        const NavSegment& s = lay.segments[i];
        if (s.start == 6 && s.end == 0 && s.kind == NavKind::Hairpin) {
            found = true;
        }
    }
    expect(found, "circular wrap segment");
}

static void test_len1_and_wrap_text() {
    NavLayout lay = must_layout(".(((..))((...)))", false);
    expect_str(lay.naview_graph, "E1(1:M0+0+0(2:P2,2:P3))", "len1 naview");
    expect_str(lay.rag_tree, "(@,@)", "len1 rag");
    expect_str(lay.rag_dual, "(@-@)", "len1 dual");

    std::string db(61, '.');
    std::string seq(61, 'N');
    NavLayout dots;
    const char* err = 0;
    expect(naview_layout(db.c_str(), 61, false, &dots, &err), "61 dots");
#if defined(_WIN32)
    _mkdir("build-linux");
    _mkdir("build-linux/naview_test_out");
#else
    mkdir("build-linux", 0755);
    mkdir("build-linux/naview_test_out", 0755);
#endif
    const char* path = "build-linux/naview_test_out/wrap.svg";
    expect(write_naview_svg(path, "wrap", seq.c_str(), db.c_str(), false, &dots, 60, 1.0, false,
                            false, 0, 0, 0),
           "wrap svg");
    expect(file_contains(path, "textLength=\"432.00\""), "full row width");
    expect(file_contains(path, "textLength=\"7.20\""), "short row width");
}

static void test_filename() {
    expect(stem_from_label("Abc.def_12", 3) == "Abc.def_12", "keep digits");
    expect(stem_from_label("a/b c", 3) == "a_b_c", "replace");
    expect(stem_from_label("...", 4) == "...", "dots kept");
    expect(stem_from_label(".", 5) == "record_5", "dot");
    expect(stem_from_label("..", 6) == "record_6", "dotdot");
    expect(stem_from_label("", 7) == "record_7", "empty");
    expect(stem_from_label("CON", 8) == "record_8", "device");
    expect(stem_from_label("com1", 9) == "record_9", "com");
    std::string long_label(100, 'A');
    long_label[90] = '/';
    std::string stem = stem_from_label(long_label, 1);
    expect(stem.size() == 80, "truncate");
    expect(stem.find('_') == std::string::npos, "slash past the cut");
}

int main() {
    test_worked_linear();
    test_worked_circular();
    test_interior_bulge();
    test_no_pairs();
    test_reject();
    test_wrap_segment();
    test_len1_and_wrap_text();
    test_filename();
    test_svg_smoke();
    test_svg_color();
    if (g_fails) {
        std::fprintf(stderr, "%d failure(s)\n", g_fails);
        return 1;
    }
    std::printf("ok\n");
    return 0;
}
