// Draw one naview SVG per row of a TSV file.
//
// rce_naview [--label-col N] [--seq-col N] [--db-col N]
//            [--out-dir DIR] [--wrap N] [--graph-width F]
//            [--rainbow | --colors LIST] [--threads N]
//            [--circular] [--header] [--naview-only] input.tsv

#include "filename.hpp"
#include "naview.hpp"
#include "svg_page.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#if defined(_WIN32)
#include <direct.h>
#include <sys/stat.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

static void die(const char* message) {
    std::fprintf(stderr, "rce_naview: %s\n", message);
    std::exit(1);
}

static void die_line(const char* path, int line_no, const char* message) {
    std::fprintf(stderr, "rce_naview: %s:%d: %s\n", path, line_no, message);
    std::exit(1);
}

static void usage(void) {
    std::fprintf(stdout,
                 "usage: rce_naview [--label-col N] [--seq-col N] [--db-col N]\n"
                 "                  [--out-dir DIR] [--wrap N] [--graph-width F]\n"
                 "                  [--rainbow | --colors LIST] [--threads N]\n"
                 "                  [--circular] [--header] [--naview-only] input.tsv\n"
                 "\n"
                 "Read label, sequence, and dot-bracket columns (1-based; default 1 2 3)\n"
                 "and write one naview SVG per record. The dot-bracket wraps at --wrap\n"
                 "symbols (default 60). The graph row is --graph-width times the naview\n"
                 "width (default 1). Groups cycle a 10-color palette. --rainbow uses one\n"
                 "hue per group, and --colors takes a comma-separated list of SVG colors.\n"
                 "--threads is the number of records laid out and drawn at once (default 1).\n"
                 "--naview-only writes the structure alone, still scaled to the\n"
                 "dot-bracket width. The SVG stem keeps [A-Za-z0-9._], replaces every\n"
                 "other byte with '_', and is truncated to 80 characters.\n");
}

static bool parse_col(const char* text, int* out) {
    if (!text || !text[0]) {
        return false;
    }
    char* end = 0;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < 1 || value > 1000000L) {
        return false;
    }
    *out = static_cast<int>(value);
    return true;
}

static bool split_colors(const char* text, std::vector<std::string>* out) {
    out->clear();
    if (!text || !text[0]) {
        return false;
    }
    std::string cur;
    int depth = 0;
    for (const char* p = text;; ++p) {
        char c = *p;
        if (c == '(') {
            depth += 1;
        } else if (c == ')' && depth > 0) {
            depth -= 1;
        }
        if ((c == ',' && depth == 0) || c == '\0') {
            size_t a = 0;
            while (a < cur.size() && (cur[a] == ' ' || cur[a] == '\t')) {
                a += 1;
            }
            size_t b = cur.size();
            while (b > a && (cur[b - 1] == ' ' || cur[b - 1] == '\t')) {
                b -= 1;
            }
            if (a == b) {
                return false;
            }
            std::string tok = cur.substr(a, b - a);
            if (!naview_svg_color(tok.c_str())) {
                return false;
            }
            out->push_back(std::move(tok));
            cur.clear();
            if (c == '\0') {
                return true;
            }
            continue;
        }
        cur.push_back(c);
    }
    return false;
}

static bool parse_factor(const char* text, double* out) {
    if (!text || !text[0]) {
        return false;
    }
    char* end = 0;
    double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || !(value > 0.0) || value > 100.0) {
        return false;
    }
    *out = value;
    return true;
}

static bool is_directory(const char* path) {
#if defined(_WIN32)
    struct _stat st;
    if (_stat(path, &st) != 0) {
        return false;
    }
    return (st.st_mode & _S_IFMT) == _S_IFDIR;
#else
    struct stat st;
    if (stat(path, &st) != 0) {
        return false;
    }
    return S_ISDIR(st.st_mode) != 0;
#endif
}

static bool path_exists(const char* path) {
#if defined(_WIN32)
    struct _stat st;
    return _stat(path, &st) == 0;
#else
    struct stat st;
    return stat(path, &st) == 0;
#endif
}

static void ensure_dir(const char* path) {
    if (path_exists(path)) {
        if (!is_directory(path)) {
            std::fprintf(stderr, "rce_naview: not a directory: %s\n", path);
            std::exit(1);
        }
        return;
    }
#if defined(_WIN32)
    if (_mkdir(path) != 0) {
#else
    if (mkdir(path, 0755) != 0) {
#endif
        std::fprintf(stderr, "rce_naview: cannot create %s: %s\n", path, std::strerror(errno));
        std::exit(1);
    }
}

static bool is_blank_or_comment(const std::string& line) {
    size_t i = 0;
    while (i < line.size() && (line[i] == ' ' || line[i] == '\t')) {
        i += 1;
    }
    return i == line.size() || line[i] == '#';
}

static bool read_line(FILE* f, std::string* out) {
    out->clear();
    int c = std::fgetc(f);
    if (c == EOF) {
        return false;
    }
    while (c != EOF && c != '\n') {
        out->push_back(static_cast<char>(c));
        c = std::fgetc(f);
    }
    if (!out->empty() && out->back() == '\r') {
        out->pop_back();
    }
    return true;
}

static std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\t') {
            fields.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(line[i]);
        }
    }
    fields.push_back(cur);
    return fields;
}

struct Record {
    std::string label;
    std::string sequence;
    std::string db;
    std::string stem;
    int line_no;
    NavLayout layout;
};

static std::string join_svg(const std::string& dir, const std::string& stem) {
    std::string path = dir;
    if (!path.empty() && path.back() != '/' && path.back() != '\\') {
        path.push_back('/');
    }
    path += stem;
    path += ".svg";
    return path;
}

struct Range {
    size_t begin;
    size_t end;
};

// One worker owns each range. A single record is never shared between workers.
static std::vector<Range> make_ranges(size_t n, int threads) {
    std::vector<Range> ranges;
    if (n == 0) {
        return ranges;
    }
    int workers = threads;
    if (static_cast<size_t>(workers) > n) {
        workers = static_cast<int>(n);
    }
    ranges.resize(static_cast<size_t>(workers));
    size_t start = 0;
    size_t base = n / static_cast<size_t>(workers);
    size_t extra = n % static_cast<size_t>(workers);
    for (int t = 0; t < workers; ++t) {
        size_t count = base + (static_cast<size_t>(t) < extra ? 1 : 0);
        ranges[static_cast<size_t>(t)].begin = start;
        ranges[static_cast<size_t>(t)].end = start + count;
        start += count;
    }
    return ranges;
}

template <typename Fn>
static void run_ranges(const std::vector<Range>& ranges, Fn fn) {
    if (ranges.empty()) {
        return;
    }
    if (ranges.size() == 1) {
        fn(0, ranges[0].begin, ranges[0].end);
        return;
    }
    std::vector<std::thread> pool;
    pool.reserve(ranges.size());
    for (size_t t = 0; t < ranges.size(); ++t) {
        pool.emplace_back(fn, static_cast<int>(t), ranges[t].begin, ranges[t].end);
    }
    for (size_t t = 0; t < pool.size(); ++t) {
        pool[t].join();
    }
}

int main(int argc, char** argv) {
    int label_col = 1;
    int seq_col = 2;
    int db_col = 3;
    const char* out_dir = ".";
    int wrap = 60;
    int threads = 1;
    double graph_width = 1.0;
    bool naview_only = false;
    bool rainbow = false;
    bool colors_set = false;
    std::vector<std::string> color_list;
    bool circular = false;
    bool header = false;
    const char* input = 0;

    for (int i = 1; i < argc; ++i) {
        const char* arg = argv[i];
        if (std::strcmp(arg, "--help") == 0) {
            usage();
            return 0;
        }
        if (std::strcmp(arg, "--circular") == 0) {
            circular = true;
            continue;
        }
        if (std::strcmp(arg, "--header") == 0) {
            header = true;
            continue;
        }
        if (std::strcmp(arg, "--naview-only") == 0) {
            naview_only = true;
            continue;
        }
        if (std::strcmp(arg, "--rainbow") == 0) {
            rainbow = true;
            continue;
        }
        const char* value = 0;
        if (std::strcmp(arg, "--label-col") == 0 || std::strcmp(arg, "--seq-col") == 0 ||
            std::strcmp(arg, "--db-col") == 0 || std::strcmp(arg, "--out-dir") == 0 ||
            std::strcmp(arg, "--wrap") == 0 || std::strcmp(arg, "--graph-width") == 0 ||
            std::strcmp(arg, "--colors") == 0 || std::strcmp(arg, "--threads") == 0) {
            if (i + 1 >= argc) {
                die("missing option value");
            }
            value = argv[++i];
        }
        if (std::strcmp(arg, "--label-col") == 0) {
            if (!parse_col(value, &label_col)) {
                die("label column must be an integer >= 1");
            }
            continue;
        }
        if (std::strcmp(arg, "--seq-col") == 0) {
            if (!parse_col(value, &seq_col)) {
                die("sequence column must be an integer >= 1");
            }
            continue;
        }
        if (std::strcmp(arg, "--db-col") == 0) {
            if (!parse_col(value, &db_col)) {
                die("dot-bracket column must be an integer >= 1");
            }
            continue;
        }
        if (std::strcmp(arg, "--out-dir") == 0) {
            if (!value || !value[0]) {
                die("output directory is empty");
            }
            out_dir = value;
            continue;
        }
        if (std::strcmp(arg, "--wrap") == 0) {
            if (!parse_col(value, &wrap)) {
                die("wrap must be an integer >= 1");
            }
            continue;
        }
        if (std::strcmp(arg, "--threads") == 0) {
            if (!parse_col(value, &threads) || threads > 256) {
                die("threads must be an integer from 1 to 256");
            }
            continue;
        }
        if (std::strcmp(arg, "--graph-width") == 0) {
            if (!parse_factor(value, &graph_width)) {
                die("graph-width must be a number > 0");
            }
            continue;
        }
        if (std::strcmp(arg, "--colors") == 0) {
            if (!split_colors(value, &color_list)) {
                die("colors must be a comma-separated list of SVG colors");
            }
            colors_set = true;
            continue;
        }
        if (arg[0] == '-' && arg[1] == '-') {
            std::fprintf(stderr, "rce_naview: unknown option %s\n", arg);
            return 1;
        }
        if (input) {
            die("only one input file is accepted");
        }
        input = arg;
    }
    if (!input) {
        usage();
        return 1;
    }
    if (label_col == seq_col || label_col == db_col || seq_col == db_col) {
        die("label, sequence, and dot-bracket columns must be distinct");
    }
    if (rainbow && colors_set) {
        die("use only one of --rainbow and --colors");
    }
    std::vector<const char*> color_ptrs;
    color_ptrs.reserve(color_list.size());
    for (size_t i = 0; i < color_list.size(); ++i) {
        color_ptrs.push_back(color_list[i].c_str());
    }

    FILE* f = std::fopen(input, "rb");
    if (!f) {
        std::fprintf(stderr, "rce_naview: cannot open %s: %s\n", input, std::strerror(errno));
        return 1;
    }
    std::vector<Record> records;
    std::string line;
    int line_no = 0;
    bool skipped_header = false;
    while (read_line(f, &line)) {
        line_no += 1;
        if (is_blank_or_comment(line)) {
            continue;
        }
        if (header && !skipped_header) {
            skipped_header = true;
            continue;
        }
        std::vector<std::string> fields = split_tabs(line);
        int need = label_col;
        if (seq_col > need) {
            need = seq_col;
        }
        if (db_col > need) {
            need = db_col;
        }
        if (static_cast<int>(fields.size()) < need) {
            std::fclose(f);
            die_line(input, line_no, "not enough columns");
        }
        Record rec;
        rec.label = fields[static_cast<size_t>(label_col - 1)];
        rec.sequence = fields[static_cast<size_t>(seq_col - 1)];
        rec.db = fields[static_cast<size_t>(db_col - 1)];
        rec.line_no = line_no;
        rec.stem = stem_from_label(rec.label, line_no);
        if (rec.sequence.empty() || rec.db.empty()) {
            std::fclose(f);
            die_line(input, line_no, "empty sequence or structure");
        }
        if (rec.sequence.size() != rec.db.size()) {
            std::fclose(f);
            die_line(input, line_no, "sequence and structure lengths differ");
        }
        records.push_back(std::move(rec));
    }
    std::fclose(f);

    struct LayoutFail {
        int line_no;
        const char* message;
    };
    std::vector<Range> ranges = make_ranges(records.size(), threads);
    std::vector<LayoutFail> layout_fails(ranges.size());
    for (size_t i = 0; i < layout_fails.size(); ++i) {
        layout_fails[i].line_no = 0;
        layout_fails[i].message = 0;
    }
    run_ranges(ranges, [&](int worker, size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            Record& rec = records[i];
            const char* err = 0;
            if (!naview_layout(rec.db.c_str(), static_cast<int>(rec.db.size()), circular,
                               &rec.layout, &err)) {
                layout_fails[static_cast<size_t>(worker)].line_no = rec.line_no;
                layout_fails[static_cast<size_t>(worker)].message = err ? err : "invalid structure";
                return;
            }
        }
    });
    int bad_line = 0;
    const char* bad_msg = 0;
    for (size_t i = 0; i < layout_fails.size(); ++i) {
        if (layout_fails[i].line_no != 0 &&
            (bad_line == 0 || layout_fails[i].line_no < bad_line)) {
            bad_line = layout_fails[i].line_no;
            bad_msg = layout_fails[i].message;
        }
    }
    if (bad_line != 0) {
        die_line(input, bad_line, bad_msg ? bad_msg : "invalid structure");
    }

    std::vector<size_t> order(records.size());
    for (size_t i = 0; i < records.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (records[a].stem != records[b].stem) {
            return records[a].stem < records[b].stem;
        }
        return records[a].line_no < records[b].line_no;
    });
    for (size_t i = 1; i < order.size(); ++i) {
        const Record& prev = records[order[i - 1]];
        const Record& cur = records[order[i]];
        if (prev.stem == cur.stem) {
            std::fprintf(stderr, "rce_naview: %s:%d and %s:%d both write %s.svg\n", input,
                         prev.line_no, input, cur.line_no, cur.stem.c_str());
            return 1;
        }
    }

    ensure_dir(out_dir);
    struct WriteFail {
        size_t index;
        int err_no;
        std::string path;
    };
    std::vector<WriteFail> write_fails(ranges.size());
    for (size_t i = 0; i < write_fails.size(); ++i) {
        write_fails[i].index = records.size();
        write_fails[i].err_no = 0;
    }
    const std::string dir(out_dir);
    const char* const* palette = color_ptrs.empty() ? 0 : color_ptrs.data();
    const int npalette = static_cast<int>(color_ptrs.size());
    run_ranges(ranges, [&](int worker, size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            const Record& rec = records[i];
            std::string path = join_svg(dir, rec.stem);
            if (!write_naview_svg(path.c_str(), rec.label.c_str(), rec.sequence.c_str(),
                                  rec.db.c_str(), circular, &rec.layout, wrap, graph_width,
                                  naview_only, rainbow, palette, npalette)) {
                write_fails[static_cast<size_t>(worker)].index = i;
                write_fails[static_cast<size_t>(worker)].err_no = errno;
                write_fails[static_cast<size_t>(worker)].path = std::move(path);
                return;
            }
        }
    });
    size_t bad_index = records.size();
    int bad_errno = 0;
    std::string bad_path;
    for (size_t i = 0; i < write_fails.size(); ++i) {
        if (write_fails[i].index < bad_index) {
            bad_index = write_fails[i].index;
            bad_errno = write_fails[i].err_no;
            bad_path = write_fails[i].path;
        }
    }
    if (bad_index < records.size()) {
        std::fprintf(stderr, "rce_naview: cannot write %s: %s\n", bad_path.c_str(),
                     std::strerror(bad_errno));
        return 1;
    }
    return 0;
}
