// Write one multi-sequence FASTA file per length/composition range.
//
// cpp_generate_fasta_grid <ranges.tsv> <output_dir> <seqs_per_range>
//
// ranges.tsv matches output/lcomp_bin_ranges.tsv: '#' comments, then
// L, fA, fC, fG, fT separated by tabs. output_dir must already exist.
// Each range is written to L<L>_<fA>_<fC>_<fG>_<fT>.fa using the TSV text.

#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>

#include <sys/stat.h>

#define asserta(cond)                                                          \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::fprintf(stderr, "asserta failed: %s:%d: %s\n", __FILE__,     \
                         __LINE__, #cond);                                     \
            std::abort();                                                      \
        }                                                                      \
    } while (0)

static const char kBases[4] = {'A', 'C', 'G', 'T'};
static const int kLineWidth = 80;
static const long long kMaxLength = 10000000LL;

static void die(const char* message) {
    std::fprintf(stderr, "cpp_generate_fasta_grid: %s\n", message);
    std::exit(1);
}

static void die_line(long long line_no, const char* message) {
    std::fprintf(stderr, "cpp_generate_fasta_grid: line %lld: %s\n", line_no, message);
    std::exit(1);
}

static int is_directory(const char* path) {
#if defined(_WIN32)
    struct _stat st;
    if (_stat(path, &st) != 0) {
        return 0;
    }
    return (st.st_mode & _S_IFMT) == _S_IFDIR;
#else
    struct stat st;
    if (stat(path, &st) != 0) {
        return 0;
    }
    return S_ISDIR(st.st_mode) != 0;
#endif
}

static int comes_first(int a, int b, const double* scaled, const int* counts, int descending) {
    double frac_a = scaled[a] - static_cast<double>(counts[a]);
    double frac_b = scaled[b] - static_cast<double>(counts[b]);
    if (frac_a != frac_b) {
        if (descending) {
            return frac_a > frac_b;
        }
        return frac_a < frac_b;
    }
    if (descending) {
        return a > b;
    }
    return a < b;
}

static void sort_order(int* order, const double* scaled, const int* counts, int descending) {
    for (int i = 1; i < 4; ++i) {
        int key = order[i];
        int j = i;
        while (j > 0 && comes_first(key, order[j - 1], scaled, counts, descending)) {
            order[j] = order[j - 1];
            --j;
        }
        order[j] = key;
    }
}

static void base_counts(long long length, const double* freqs, int* counts) {
    double total = freqs[0] + freqs[1] + freqs[2] + freqs[3];
    asserta(total > 0.0);
    double scaled[4];
    int floored = 0;
    for (int i = 0; i < 4; ++i) {
        scaled[i] = freqs[i] / total * static_cast<double>(length);
        counts[i] = static_cast<int>(std::floor(scaled[i]));
        asserta(counts[i] >= 0);
        floored += counts[i];
    }
    int leftover = static_cast<int>(length) - floored;
    int order[4] = {0, 1, 2, 3};
    if (leftover >= 0) {
        sort_order(order, scaled, counts, 1);
        asserta(leftover <= 4);
        for (int i = 0; i < leftover; ++i) {
            counts[order[i]] += 1;
        }
    } else {
        sort_order(order, scaled, counts, 0);
        for (int i = 0; i < 4; ++i) {
            if (leftover == 0) {
                break;
            }
            int index = order[i];
            if (counts[index] == 0) {
                continue;
            }
            counts[index] -= 1;
            leftover += 1;
        }
        asserta(leftover == 0);
    }
    long long sum = 0;
    for (int i = 0; i < 4; ++i) {
        asserta(counts[i] >= 0);
        asserta(std::fabs(static_cast<double>(counts[i]) - scaled[i]) < 1.0);
        sum += counts[i];
    }
    asserta(sum == length);
}

// Unbiased integer in [0, bound). bound is at least 1 and fits in 32 bits.
static std::uint32_t uniform_below(std::mt19937& rng, std::uint32_t bound) {
    std::uint32_t x = rng();
    std::uint64_t product = static_cast<std::uint64_t>(x) * static_cast<std::uint64_t>(bound);
    std::uint32_t low = static_cast<std::uint32_t>(product);
    if (low < bound) {
        std::uint32_t threshold = static_cast<std::uint32_t>(-bound) % bound;
        while (low < threshold) {
            x = rng();
            product = static_cast<std::uint64_t>(x) * static_cast<std::uint64_t>(bound);
            low = static_cast<std::uint32_t>(product);
        }
    }
    return static_cast<std::uint32_t>(product >> 32);
}

static void fill_sequence(char* seq, long long length, const int* counts) {
    size_t pos = 0;
    for (int base = 0; base < 4; ++base) {
        for (int k = 0; k < counts[base]; ++k) {
            seq[pos] = kBases[base];
            pos += 1;
        }
    }
    asserta(pos == static_cast<size_t>(length));
}

static void shuffle_sequence(char* seq, long long length, long long seed) {
    asserta(seed >= 0);
    asserta(seed <= static_cast<long long>(UINT32_MAX));
    std::mt19937 rng(static_cast<std::uint32_t>(seed));
    for (long long i = length - 1; i > 0; --i) {
        std::uint32_t j = uniform_below(rng, static_cast<std::uint32_t>(i + 1));
        char tmp = seq[i];
        seq[i] = seq[j];
        seq[j] = tmp;
    }
}

static int write_record(std::FILE* out, const char* range_name, long long seed, const char* seq, long long length) {
    if (std::fprintf(out, ">%s_seed%lld\n", range_name, seed) < 0) {
        return 0;
    }
    long long offset = 0;
    while (offset < length) {
        long long chunk = length - offset;
        if (chunk > kLineWidth) {
            chunk = kLineWidth;
        }
        size_t n = std::fwrite(seq + offset, 1, static_cast<size_t>(chunk), out);
        if (n != static_cast<size_t>(chunk) || std::fputc('\n', out) != '\n') {
            return 0;
        }
        offset += chunk;
    }
    return 1;
}

static int field_ok(const char* field) {
    if (field[0] == '\0') {
        return 0;
    }
    for (const char* p = field; *p != '\0'; ++p) {
        if (*p == '/' || *p == '\\' || *p == ':') {
            return 0;
        }
    }
    return 1;
}

static long long parse_i64(const char* text, int* ok) {
    char* end = nullptr;
    errno = 0;
    long long value = std::strtoll(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0') {
        *ok = 0;
        return 0;
    }
    *ok = 1;
    return value;
}

static double parse_f64(const char* text, int* ok) {
    char* end = nullptr;
    errno = 0;
    double value = std::strtod(text, &end);
    if (errno != 0 || end == text || *end != '\0' || !std::isfinite(value)) {
        *ok = 0;
        return 0.0;
    }
    *ok = 1;
    return value;
}

static void process_range(const char* out_dir, const char* fields[5], long long length, const double* freqs, long long n_seqs, long long line_no) {
    for (int i = 0; i < 5; ++i) {
        if (!field_ok(fields[i])) {
            die_line(line_no, "range field is empty or not a single path component");
        }
    }
    char range_name[256];
    int name_n = std::snprintf(range_name, sizeof range_name, "L%s_%s_%s_%s_%s", fields[0], fields[1], fields[2], fields[3], fields[4]);
    if (name_n < 0 || static_cast<size_t>(name_n) >= sizeof range_name) {
        die("range name is too long");
    }

    int need_slash = 1;
    size_t dir_len = std::strlen(out_dir);
    if (dir_len > 0 && (out_dir[dir_len - 1] == '/' || out_dir[dir_len - 1] == '\\')) {
        need_slash = 0;
    }
    char path[4096];
    int path_n = std::snprintf(path, sizeof path, need_slash ? "%s/%s.fa" : "%s%s.fa", out_dir, range_name);
    if (path_n < 0 || static_cast<size_t>(path_n) >= sizeof path) {
        die("output path is too long");
    }

    std::FILE* out = std::fopen(path, "wb");
    if (out == nullptr) {
        std::fprintf(stderr, "cpp_generate_fasta_grid: cannot write %s\n", path);
        std::exit(1);
    }
    char* io_buf = static_cast<char*>(std::malloc(1 << 20));
    if (io_buf != nullptr) {
        std::setvbuf(out, io_buf, _IOFBF, 1 << 20);
    }

    int counts[4];
    base_counts(length, freqs, counts);
    char* seq = nullptr;
    if (n_seqs > 0) {
        asserta(length > 0);
        seq = static_cast<char*>(std::malloc(static_cast<size_t>(length)));
        if (seq == nullptr) {
            die("out of memory");
        }
    }
    for (long long seed = 0; seed < n_seqs; ++seed) {
        fill_sequence(seq, length, counts);
        shuffle_sequence(seq, length, seed);
        if (!write_record(out, range_name, seed, seq, length) || std::ferror(out)) {
            std::fprintf(stderr, "cpp_generate_fasta_grid: write failed for %s\n", path);
            std::exit(1);
        }
    }
    if (std::fclose(out) != 0) {
        std::fprintf(stderr, "cpp_generate_fasta_grid: close failed for %s\n", path);
        std::exit(1);
    }
    std::free(seq);
    std::free(io_buf);
    std::printf("%s\n", path);
}

static void read_ranges(const char* tsv_path, const char* out_dir, long long n_seqs) {
    std::FILE* in = std::fopen(tsv_path, "rb");
    if (in == nullptr) {
        std::fprintf(stderr, "cpp_generate_fasta_grid: cannot read %s\n", tsv_path);
        std::exit(1);
    }
    char line[4096];
    long long line_no = 0;
    long long n_files = 0;
    while (std::fgets(line, static_cast<int>(sizeof line), in) != nullptr) {
        line_no += 1;
        size_t len = std::strlen(line);
        if (len == 0) {
            continue;
        }
        if (line[len - 1] != '\n' && !std::feof(in)) {
            die_line(line_no, "line is too long");
        }
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r')) {
            line[len - 1] = '\0';
            len -= 1;
        }
        if (len == 0 || line[0] == '#') {
            continue;
        }

        const char* fields[5];
        int n_fields = 0;
        int extra = 0;
        char* cursor = line;
        while (n_fields < 5) {
            fields[n_fields] = cursor;
            n_fields += 1;
            char* tab = std::strchr(cursor, '\t');
            if (tab == nullptr) {
                break;
            }
            *tab = '\0';
            cursor = tab + 1;
            if (n_fields == 5) {
                extra = 1;
            }
        }
        if (n_fields != 5 || extra) {
            die_line(line_no, "expected 5 tab-separated fields");
        }

        int ok = 0;
        long long length = parse_i64(fields[0], &ok);
        if (!ok || length < 1 || length > kMaxLength) {
            die_line(line_no, "L must be an integer from 1 to 10000000");
        }
        double freqs[4];
        for (int i = 0; i < 4; ++i) {
            freqs[i] = parse_f64(fields[i + 1], &ok);
            if (!ok || freqs[i] < 0.0) {
                die_line(line_no, "frequency must be a finite number >= 0");
            }
        }
        double total = freqs[0] + freqs[1] + freqs[2] + freqs[3];
        if (total < 0.99 || total > 1.01) {
            die_line(line_no, "frequencies sum outside 0.99..1.01");
        }
        process_range(out_dir, fields, length, freqs, n_seqs, line_no);
        n_files += 1;
    }
    if (std::ferror(in)) {
        die("error reading ranges");
    }
    std::fclose(in);
    std::fprintf(stderr, "wrote %lld files, %lld sequences each\n", n_files, n_seqs);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::fprintf(stderr, "usage: cpp_generate_fasta_grid <ranges.tsv> <output_dir> <seqs_per_range>\n");
        return 1;
    }
    const char* tsv_path = argv[1];
    const char* out_dir = argv[2];
    int ok = 0;
    long long n_seqs = parse_i64(argv[3], &ok);
    if (!ok || n_seqs < 0) {
        die("seqs_per_range must be an integer >= 0");
    }
    if (n_seqs > static_cast<long long>(UINT32_MAX) + 1LL) {
        die("seqs_per_range is too large");
    }
    if (!is_directory(out_dir)) {
        std::fprintf(stderr, "cpp_generate_fasta_grid: output directory does not exist: %s\n", out_dir);
        return 1;
    }
    read_ranges(tsv_path, out_dir, n_seqs);
    return 0;
}
