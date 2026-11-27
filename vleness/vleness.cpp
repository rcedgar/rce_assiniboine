// Score circular RNA dot-brackets with the published VLE formula
// S = (pct_bulge - 12) / 5.21.
//
// vleness [--threads N] [--rnafold | --tsv I,J | --lsc] [--keep-full-label]
//         [--log FILE] [--maxerror N] [--minscore X] input output.tsv
//
// input may be a path, including /dev/stdin. The file is never seeked.

#include "vleness.h"

#include "../rce_naview/asserta.hpp"

#include <cctype>
#include <clocale>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <thread>
#include <vector>

static const int kGap = 5;
static const double kBulgeMedian = 12.0;
static const double kBulgeIqr = 5.21;
static const int kDefaultMaxError = 1000;
static const double kDefaultMinScore = 0.75;

static const char* kUsage =
    "usage: vleness [--threads N] [--rnafold | --tsv I,J | --lsc]\n"
    "               [--keep-full-label] [--log FILE] [--maxerror N]\n"
    "               [--minscore X] input output.tsv\n"
    "\n"
    "Score each circular dot-bracket. pct_bulge is the percent of the circle\n"
    "in unpaired gaps inside composite stems (helices joined across at most\n"
    "five unpaired nucleotides on a strand). S = (pct_bulge - 12) / 5.21.\n"
    "pct_bulge is printed with one decimal place and S with four.\n"
    "Rows with S below --minscore are omitted (default 0.75).\n"
    "An LSC row is kept when the larger of S1 and S2 reaches that cutoff.\n"
    "\n"
    "The first byte selects the format: '>' is RNAfold, anything else is TSV.\n"
    "--rnafold forces three-line RNAfold records (label, sequence, structure).\n"
    "A trailing Vienna energy on the structure line is ignored. The label is\n"
    "cut at the first whitespace unless --keep-full-label is set.\n"
    "--tsv I,J reads a TSV whose 1-based columns are the label and the\n"
    "dot-bracket (default 1,2). A first row whose structure field is not a\n"
    "dot-bracket is a header and is skipped.\n"
    "--lsc reads label, MFE, dot-bracket, MFE, dot-bracket and writes\n"
    "label, length, pct_bulge1, pct_bulge2, S1, S2, max(S1, S2) with no header.\n"
    "LSC is also chosen automatically when the first TSV row has dot-brackets\n"
    "in fields 3 and 5. --rnafold, --tsv, and --lsc are mutually exclusive.\n"
    "\n"
    "Format problems are warnings on stderr. --log FILE also records them,\n"
    "with the start time and the finish time. The run stops after more than\n"
    "--maxerror warnings (default 1000) and exits non-zero. Rows accepted\n"
    "before that stop are still written.\n"
    "--threads defaults to the processor count in /proc/cpuinfo.\n";

enum class Format {
    Auto,
    Rnafold,
    Tsv,
    Lsc,
};

struct Options {
    Format format;
    int label_col;
    int db_col;
    int threads;
    int max_error;
    double min_score;
    bool keep_full_label;
    const char* log_path;
    const char* input_path;
    const char* output_path;
};

struct Job {
    std::string label;
    std::string db1;
    std::string db2;
    int length;
    double pct1;
    double pct2;
    double s1;
    double s2;
    bool lsc;
};

static FILE* g_log = 0;

static void format_now(char* buf, size_t n) {
    std::time_t now = std::time(0);
    std::tm local;
#if defined(_WIN32)
    localtime_s(&local, &now);
#else
    localtime_r(&now, &local);
#endif
    if (std::strftime(buf, n, "%Y-%m-%d %H:%M:%S", &local) == 0) {
        buf[0] = '\0';
    }
}

static void write_finish(void) {
    if (!g_log) {
        return;
    }
    char when[32];
    format_now(when, sizeof when);
    std::fprintf(g_log, "finish\t%s\n", when);
    std::fclose(g_log);
    g_log = 0;
}

static void usage(FILE* out) {
    std::fputs(kUsage, out);
}

static void die_usage(const char* message) {
    std::fprintf(stderr, "vleness: %s\n", message);
    usage(stderr);
    std::exit(1);
}

static bool parse_i64(const char* text, long* out, long lo, long hi) {
    if (!text || !text[0]) {
        return false;
    }
    char* end = 0;
    long value = std::strtol(text, &end, 10);
    if (end == text || *end != '\0' || value < lo || value > hi) {
        return false;
    }
    *out = value;
    return true;
}

static bool parse_f64(const char* text, double* out) {
    if (!text || !text[0]) {
        return false;
    }
    char* end = 0;
    double value = std::strtod(text, &end);
    if (end == text || *end != '\0' || !std::isfinite(value)) {
        return false;
    }
    *out = value;
    return true;
}

static bool parse_tsv_spec(const char* text, int* label_col, int* db_col) {
    if (!text || !text[0]) {
        return false;
    }
    const char* comma = std::strchr(text, ',');
    if (!comma || std::strchr(comma + 1, ',')) {
        return false;
    }
    std::string left(text, comma);
    long a = 0;
    long b = 0;
    if (!parse_i64(left.c_str(), &a, 1, 1000000L)) {
        return false;
    }
    if (!parse_i64(comma + 1, &b, 1, 1000000L)) {
        return false;
    }
    *label_col = static_cast<int>(a);
    *db_col = static_cast<int>(b);
    return true;
}

static int default_threads(void) {
    FILE* cpuinfo = std::fopen("/proc/cpuinfo", "r");
    if (cpuinfo) {
        int count = 0;
        char line[256];
        while (std::fgets(line, sizeof line, cpuinfo)) {
            if (std::strncmp(line, "processor", 9) == 0) {
                count += 1;
            }
        }
        std::fclose(cpuinfo);
        if (count > 0) {
            return count;
        }
    }
    unsigned hardware = std::thread::hardware_concurrency();
    if (hardware == 0) {
        return 1;
    }
    return static_cast<int>(hardware);
}

static bool read_line(FILE* in, std::string* line) {
    line->clear();
    int ch = std::fgetc(in);
    if (ch == EOF) {
        return false;
    }
    while (ch != EOF && ch != '\n') {
        line->push_back(static_cast<char>(ch));
        ch = std::fgetc(in);
    }
    return true;
}

static int modn(int value, int n) {
    int r = value % n;
    if (r < 0) {
        r += n;
    }
    return r;
}

static bool parse_pairs(const std::string& db, std::vector<int>* partner, const char** error) {
    const int n = static_cast<int>(db.size());
    if (n <= 0) {
        *error = "empty dot-bracket";
        return false;
    }
    partner->assign(static_cast<size_t>(n), -1);
    std::vector<int> stack;
    stack.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        char c = db[static_cast<size_t>(i)];
        if (c == '(') {
            stack.push_back(i);
        } else if (c == ')') {
            if (stack.empty()) {
                *error = "unmatched closing bracket";
                return false;
            }
            int j = stack.back();
            stack.pop_back();
            (*partner)[static_cast<size_t>(j)] = i;
            (*partner)[static_cast<size_t>(i)] = j;
        } else if (c != '.') {
            *error = "invalid character";
            return false;
        }
    }
    if (!stack.empty()) {
        *error = "unmatched opening bracket";
        return false;
    }
    *error = 0;
    return true;
}

static bool next_inward(
    const std::string& db,
    const std::vector<int>& partner,
    int i,
    int j,
    int* ni,
    int* nj,
    int* len5,
    int* len3
) {
    const int n = static_cast<int>(db.size());
    for (int step5 = 1; step5 <= kGap + 1; ++step5) {
        int cand_i = modn(i + step5, n);
        if (db[static_cast<size_t>(cand_i)] == '.') {
            continue;
        }
        for (int step3 = 1; step3 <= kGap + 1; ++step3) {
            int cand_j = modn(j - step3, n);
            if (db[static_cast<size_t>(cand_j)] == '.') {
                continue;
            }
            if (partner[static_cast<size_t>(cand_i)] == cand_j && cand_i != j && cand_j != i) {
                *ni = cand_i;
                *nj = cand_j;
                *len5 = step5 - 1;
                *len3 = step3 - 1;
                return true;
            }
            return false;
        }
        return false;
    }
    return false;
}

static bool prev_outward(
    const std::string& db,
    const std::vector<int>& partner,
    int i,
    int j,
    int* pi,
    int* pj
) {
    const int n = static_cast<int>(db.size());
    for (int step5 = 1; step5 <= kGap + 1; ++step5) {
        int cand_i = modn(i - step5, n);
        if (db[static_cast<size_t>(cand_i)] == '.') {
            continue;
        }
        for (int step3 = 1; step3 <= kGap + 1; ++step3) {
            int cand_j = modn(j + step3, n);
            if (db[static_cast<size_t>(cand_j)] == '.') {
                continue;
            }
            if (partner[static_cast<size_t>(cand_i)] == cand_j && cand_i != j && cand_j != i) {
                *pi = cand_i;
                *pj = cand_j;
                return true;
            }
            return false;
        }
        return false;
    }
    return false;
}

static int bulge_bases(const std::string& db, const std::vector<int>& partner) {
    const int n = static_cast<int>(db.size());
    std::vector<char> visited(static_cast<size_t>(n), 0);
    int total = 0;
    for (int i = 0; i < n; ++i) {
        int j = partner[static_cast<size_t>(i)];
        if (j <= i) {
            continue;
        }
        if (visited[static_cast<size_t>(i)] || visited[static_cast<size_t>(j)]) {
            continue;
        }
        int ci = i;
        int cj = j;
        while (true) {
            int pi = 0;
            int pj = 0;
            if (!prev_outward(db, partner, ci, cj, &pi, &pj)) {
                break;
            }
            if (visited[static_cast<size_t>(pi)] || visited[static_cast<size_t>(pj)]) {
                break;
            }
            if (pi == i && pj == j) {
                break;
            }
            ci = pi;
            cj = pj;
        }
        while (!visited[static_cast<size_t>(ci)] && !visited[static_cast<size_t>(cj)]) {
            visited[static_cast<size_t>(ci)] = 1;
            visited[static_cast<size_t>(cj)] = 1;
            int ni = 0;
            int nj = 0;
            int len5 = 0;
            int len3 = 0;
            if (!next_inward(db, partner, ci, cj, &ni, &nj, &len5, &len3)) {
                break;
            }
            if (len5 > 0) {
                total += len5;
            }
            if (len3 > 0) {
                total += len3;
            }
            ci = ni;
            cj = nj;
        }
    }
    return total;
}

bool vleness_score(const std::string& dot_bracket, VlenessScore* score, std::string* error) {
    asserta(score != 0);
    const char* problem = 0;
    std::vector<int> partner;
    if (!parse_pairs(dot_bracket, &partner, &problem)) {
        if (error) {
            *error = problem;
        }
        return false;
    }
    const int n = static_cast<int>(dot_bracket.size());
    const int bases = bulge_bases(dot_bracket, partner);
    score->length = n;
    score->pct_bulge = 100.0 * static_cast<double>(bases) / static_cast<double>(n);
    score->score = (score->pct_bulge - kBulgeMedian) / kBulgeIqr;
    return true;
}

static bool is_dot_bracket(const std::string& text) {
    if (text.empty()) {
        return false;
    }
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        if (c != '.' && c != '(' && c != ')') {
            return false;
        }
    }
    return true;
}

static void strip_cr(std::string* line) {
    if (!line->empty() && line->back() == '\r') {
        line->pop_back();
    }
}

static void split_tabs(const std::string& line, std::vector<std::string>* fields) {
    fields->clear();
    std::string cur;
    for (size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '\t') {
            fields->push_back(cur);
            cur.clear();
        } else {
            cur.push_back(line[i]);
        }
    }
    fields->push_back(cur);
}

static std::string field_at(const std::vector<std::string>& fields, int col_1based) {
    if (col_1based < 1 || static_cast<size_t>(col_1based) > fields.size()) {
        return std::string();
    }
    return fields[static_cast<size_t>(col_1based - 1)];
}

static bool has_field(const std::vector<std::string>& fields, int col_1based) {
    return col_1based >= 1 && static_cast<size_t>(col_1based) <= fields.size();
}

static std::string rnafold_label(const std::string& header, bool keep_full) {
    std::string label = header;
    if (!label.empty() && label[0] == '>') {
        label.erase(label.begin());
    }
    if (!keep_full) {
        size_t cut = 0;
        while (cut < label.size() && !std::isspace(static_cast<unsigned char>(label[cut]))) {
            cut += 1;
        }
        label.resize(cut);
    }
    return label;
}

static std::string structure_token(const std::string& line) {
    size_t end = 0;
    while (end < line.size() && !std::isspace(static_cast<unsigned char>(line[end]))) {
        end += 1;
    }
    return line.substr(0, end);
}

struct Reader {
    Options opt;
    FILE* in;
    std::vector<Job>* jobs;
    int errors;
    bool stop;
    bool lsc;
    long line_no;

    bool warn(const std::string& label, long line, const char* message) {
        const char* shown = label.empty() ? "-" : label.c_str();
        std::fprintf(stderr, "warning\t%s\tline %ld: %s\n", shown, line, message);
        if (g_log) {
            std::fprintf(g_log, "warning\t%s\tline %ld: %s\n", shown, line, message);
            std::fflush(g_log);
        }
        errors += 1;
        if (errors > opt.max_error) {
            std::fprintf(
                stderr, "fatal\terror count %d exceeds maxerror %d\n", errors, opt.max_error
            );
            if (g_log) {
                std::fprintf(
                    g_log, "fatal\terror count %d exceeds maxerror %d\n", errors, opt.max_error
                );
                std::fflush(g_log);
            }
            stop = true;
            return false;
        }
        return true;
    }

    bool accept(const std::string& label, const std::string& db1, const std::string& db2, bool lsc_row, long line) {
        const char* problem = 0;
        std::vector<int> partner;
        if (!parse_pairs(db1, &partner, &problem)) {
            return warn(label, line, problem);
        }
        if (lsc_row) {
            if (db2.size() != db1.size()) {
                return warn(label, line, "structure length differs");
            }
            if (!parse_pairs(db2, &partner, &problem)) {
                return warn(label, line, problem);
            }
        }
        Job job;
        job.label = label;
        job.db1 = db1;
        job.db2 = db2;
        job.length = static_cast<int>(db1.size());
        job.pct1 = 0.0;
        job.s1 = 0.0;
        job.pct2 = 0.0;
        job.s2 = 0.0;
        job.lsc = lsc_row;
        jobs->push_back(job);
        return true;
    }
};

static bool read_rnafold(Reader* reader, const std::string* first_line) {
    std::string header;
    std::string sequence;
    std::string structure;
    bool have_first = first_line != 0;
    while (!reader->stop) {
        if (have_first) {
            header = *first_line;
            have_first = false;
        } else if (!read_line(reader->in, &header)) {
            break;
        }
        reader->line_no += 1;
        strip_cr(&header);
        long header_line = reader->line_no;
        if (header.empty()) {
            if (!reader->warn("", header_line, "empty line")) {
                return false;
            }
            continue;
        }
        if (header[0] != '>') {
            if (!reader->warn("", header_line, "expected a FASTA label")) {
                return false;
            }
            continue;
        }
        std::string label = rnafold_label(header, reader->opt.keep_full_label);
        if (!read_line(reader->in, &sequence)) {
            reader->warn(label, header_line, "truncated record");
            return !reader->stop;
        }
        reader->line_no += 1;
        strip_cr(&sequence);
        if (!read_line(reader->in, &structure)) {
            reader->warn(label, header_line, "truncated record");
            return !reader->stop;
        }
        reader->line_no += 1;
        strip_cr(&structure);
        std::string db = structure_token(structure);
        if (db.size() != sequence.size()) {
            if (!reader->warn(label, header_line, "sequence length differs")) {
                return false;
            }
            continue;
        }
        if (!reader->accept(label, db, "", false, header_line)) {
            return false;
        }
    }
    return !reader->stop;
}

static bool looks_like_lsc(const std::vector<std::string>& fields) {
    if (fields.size() < 5) {
        return false;
    }
    return is_dot_bracket(fields[2]) && is_dot_bracket(fields[4]);
}

static bool read_tsv(Reader* reader, std::string first_line) {
    strip_cr(&first_line);
    std::vector<std::string> fields;
    split_tabs(first_line, &fields);
    reader->line_no = 1;
    bool explicit_tsv = reader->opt.format == Format::Tsv;
    if (reader->opt.format == Format::Lsc) {
        reader->lsc = true;
    } else if (!explicit_tsv && looks_like_lsc(fields)) {
        reader->lsc = true;
    }
    int label_col = reader->lsc ? 1 : reader->opt.label_col;
    int db_col = reader->lsc ? 3 : reader->opt.db_col;
    int db2_col = 5;
    bool header = false;
    if (reader->lsc) {
        if (has_field(fields, 3) && has_field(fields, 5) &&
            (!is_dot_bracket(field_at(fields, 3)) || !is_dot_bracket(field_at(fields, 5)))) {
            header = true;
        }
    } else if (has_field(fields, db_col) && !is_dot_bracket(field_at(fields, db_col))) {
        header = true;
    }
    bool use_first = !header;
    while (!reader->stop) {
        if (!use_first) {
            if (!read_line(reader->in, &first_line)) {
                break;
            }
            reader->line_no += 1;
            strip_cr(&first_line);
            split_tabs(first_line, &fields);
        }
        use_first = false;
        long line = reader->line_no;
        if (first_line.empty()) {
            if (!reader->warn("", line, "empty line")) {
                return false;
            }
            continue;
        }
        std::string label;
        if (!has_field(fields, label_col) || !has_field(fields, db_col) ||
            (reader->lsc && !has_field(fields, db2_col))) {
            if (!reader->warn(field_at(fields, label_col), line, "missing column")) {
                return false;
            }
            continue;
        }
        label = field_at(fields, label_col);
        std::string db1 = field_at(fields, db_col);
        std::string db2 = reader->lsc ? field_at(fields, db2_col) : std::string();
        if (!reader->accept(label, db1, db2, reader->lsc, line)) {
            return false;
        }
    }
    return !reader->stop;
}

static void score_range(std::vector<Job>* jobs, size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
        Job* job = &(*jobs)[i];
        VlenessScore score;
        std::string error;
        asserta(vleness_score(job->db1, &score, &error));
        job->length = score.length;
        job->pct1 = score.pct_bulge;
        job->s1 = score.score;
        if (job->lsc) {
            asserta(vleness_score(job->db2, &score, &error));
            job->pct2 = score.pct_bulge;
            job->s2 = score.score;
        }
    }
}

static bool keep_job(const Job& job, double min_score) {
    double best = job.s1;
    if (job.lsc && job.s2 > best) {
        best = job.s2;
    }
    return best >= min_score;
}

static bool write_output(const char* path, const std::vector<Job>& jobs, bool lsc, double min_score) {
    FILE* out = std::fopen(path, "wb");
    if (!out) {
        std::fprintf(stderr, "vleness: cannot write %s\n", path);
        return false;
    }
    if (!lsc) {
        std::fputs("label\tlength\tpct_bulge\tS\n", out);
    }
    for (size_t i = 0; i < jobs.size(); ++i) {
        const Job& job = jobs[i];
        if (!keep_job(job, min_score)) {
            continue;
        }
        if (job.lsc) {
            double best = job.s1 > job.s2 ? job.s1 : job.s2;
            std::fprintf(
                out,
                "%s\t%d\t%.1f\t%.1f\t%.4f\t%.4f\t%.4f\n",
                job.label.c_str(),
                job.length,
                job.pct1,
                job.pct2,
                job.s1,
                job.s2,
                best
            );
        } else {
            std::fprintf(
                out,
                "%s\t%d\t%.1f\t%.4f\n",
                job.label.c_str(),
                job.length,
                job.pct1,
                job.s1
            );
        }
    }
    if (std::fclose(out) != 0) {
        std::fprintf(stderr, "vleness: cannot finish %s\n", path);
        return false;
    }
    return true;
}

static void parse_args(int argc, char** argv, Options* opt) {
    opt->format = Format::Auto;
    opt->label_col = 1;
    opt->db_col = 2;
    opt->threads = default_threads();
    opt->max_error = kDefaultMaxError;
    opt->min_score = kDefaultMinScore;
    opt->keep_full_label = false;
    opt->log_path = 0;
    opt->input_path = 0;
    opt->output_path = 0;
    int positionals = 0;
    const char* positional[2] = {0, 0};
    for (int i = 1; i < argc; ++i) {
        const char* arg = argv[i];
        if (std::strcmp(arg, "-h") == 0 || std::strcmp(arg, "--help") == 0) {
            usage(stdout);
            std::exit(0);
        }
        if (std::strcmp(arg, "--rnafold") == 0) {
            if (opt->format != Format::Auto) {
                die_usage("--rnafold, --tsv, and --lsc are mutually exclusive");
            }
            opt->format = Format::Rnafold;
            continue;
        }
        if (std::strcmp(arg, "--lsc") == 0) {
            if (opt->format != Format::Auto) {
                die_usage("--rnafold, --tsv, and --lsc are mutually exclusive");
            }
            opt->format = Format::Lsc;
            continue;
        }
        if (std::strcmp(arg, "--keep-full-label") == 0) {
            opt->keep_full_label = true;
            continue;
        }
        const char* value = 0;
        const char* name = arg;
        if (std::strncmp(arg, "--tsv", 5) == 0 || std::strncmp(arg, "--threads", 9) == 0 ||
            std::strncmp(arg, "--log", 5) == 0 || std::strncmp(arg, "--maxerror", 10) == 0 ||
            std::strncmp(arg, "--minscore", 10) == 0) {
            const char* eq = std::strchr(arg, '=');
            if (eq) {
                value = eq + 1;
            } else if (i + 1 < argc) {
                value = argv[++i];
            } else {
                die_usage("missing option value");
            }
            if (std::strncmp(arg, "--tsv", 5) == 0) {
                name = "--tsv";
            } else if (std::strncmp(arg, "--threads", 9) == 0) {
                name = "--threads";
            } else if (std::strncmp(arg, "--maxerror", 10) == 0) {
                name = "--maxerror";
            } else if (std::strncmp(arg, "--minscore", 10) == 0) {
                name = "--minscore";
            } else {
                name = "--log";
            }
        }
        if (name != arg || value) {
            if (std::strcmp(name, "--tsv") == 0) {
                if (opt->format != Format::Auto) {
                    die_usage("--rnafold, --tsv, and --lsc are mutually exclusive");
                }
                if (!parse_tsv_spec(value, &opt->label_col, &opt->db_col)) {
                    die_usage("--tsv needs two 1-based columns, as 1,2");
                }
                opt->format = Format::Tsv;
                continue;
            }
            if (std::strcmp(name, "--threads") == 0) {
                long n = 0;
                if (!parse_i64(value, &n, 1, 1000000L)) {
                    die_usage("--threads must be an integer >= 1");
                }
                opt->threads = static_cast<int>(n);
                continue;
            }
            if (std::strcmp(name, "--maxerror") == 0) {
                long n = 0;
                if (!parse_i64(value, &n, 0, 1000000000L)) {
                    die_usage("--maxerror must be an integer >= 0");
                }
                opt->max_error = static_cast<int>(n);
                continue;
            }
            if (std::strcmp(name, "--minscore") == 0) {
                double score = 0.0;
                if (!parse_f64(value, &score)) {
                    die_usage("--minscore must be a number");
                }
                opt->min_score = score;
                continue;
            }
            if (std::strcmp(name, "--log") == 0) {
                if (!value || !value[0]) {
                    die_usage("--log needs a file name");
                }
                opt->log_path = value;
                continue;
            }
        }
        if (arg[0] == '-') {
            die_usage("unknown option");
        }
        if (positionals >= 2) {
            die_usage("too many file arguments");
        }
        positional[positionals++] = arg;
    }
    if (positionals != 2) {
        die_usage("need an input file and an output file");
    }
    opt->input_path = positional[0];
    opt->output_path = positional[1];
}

int main(int argc, char** argv) {
    std::setlocale(LC_ALL, "C");
    if (argc <= 1) {
        usage(stdout);
        return 0;
    }
    Options opt;
    parse_args(argc, argv, &opt);
    if (opt.log_path) {
        g_log = std::fopen(opt.log_path, "wb");
        if (!g_log) {
            std::fprintf(stderr, "vleness: cannot write %s\n", opt.log_path);
            return 1;
        }
        char when[32];
        format_now(when, sizeof when);
        std::fprintf(g_log, "start\t%s\n", when);
        std::fflush(g_log);
        std::atexit(write_finish);
    }
    FILE* input = std::fopen(opt.input_path, "rb");
    if (!input) {
        std::fprintf(stderr, "vleness: cannot read %s\n", opt.input_path);
        return 1;
    }
    std::string first;
    bool any = read_line(input, &first);
    if (any) {
        strip_cr(&first);
    }
    bool rnafold = opt.format == Format::Rnafold;
    if (opt.format == Format::Auto && any && !first.empty() && first[0] == '>') {
        rnafold = true;
    }
    std::vector<Job> jobs;
    Reader reader;
    reader.opt = opt;
    reader.in = input;
    reader.jobs = &jobs;
    reader.errors = 0;
    reader.stop = false;
    reader.lsc = opt.format == Format::Lsc;
    reader.line_no = 0;
    bool ok = true;
    if (any) {
        if (rnafold) {
            ok = read_rnafold(&reader, &first);
        } else {
            ok = read_tsv(&reader, first);
        }
    }
    int workers = opt.threads;
    if (workers < 1) {
        workers = 1;
    }
    if (!jobs.empty()) {
        if (workers > static_cast<int>(jobs.size())) {
            workers = static_cast<int>(jobs.size());
        }
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(workers));
        size_t begin = 0;
        for (int t = 0; t < workers; ++t) {
            size_t end = jobs.size() * static_cast<size_t>(t + 1) / static_cast<size_t>(workers);
            pool.emplace_back(score_range, &jobs, begin, end);
            begin = end;
        }
        for (size_t t = 0; t < pool.size(); ++t) {
            pool[t].join();
        }
    }
    std::fclose(input);
    if (!write_output(opt.output_path, jobs, reader.lsc, opt.min_score)) {
        return 1;
    }
    if (!ok) {
        return 1;
    }
    return 0;
}
