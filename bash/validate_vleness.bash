#!/usr/bin/env bash
# Build vleness and compare it with extract_circular_rna_features on the
# 100-per-family label subset. Also round-trip RNAfold and LSC fixtures.
set -euo pipefail
cd "$(dirname "$0")/.."
export LC_ALL=C
mkdir -p validate_vleness

make -C vleness
bin=vleness/vleness

"$bin" > validate_vleness/usage.txt
"$bin" -h > /dev/null
"$bin" --help > /dev/null
if "$bin" --not-an-option in out > /dev/null 2> validate_vleness/bad_usage.txt; then
    echo "invalid command line exited 0" >&2
    exit 1
fi

python3 - <<'PY'
import sys
sys.path.insert(0, "py")
from dotbracket_circular_features import extract_circular_rna_features

rows = [
    ("hairpin", "(((...)))"),
    ("bulge", "(.(...))"),
    ("internal", "(.(.).)"),
    ("dots", "......"),
    ("stack", "(((())))"),
]
with open("validate_vleness/quick.tsv", "w", newline="\n") as raw, \
     open("validate_vleness/quick_expected.tsv", "w", newline="\n") as exp:
    exp.write("label\tlength\tpct_bulge\tS\n")
    for label, db in rows:
        raw.write("%s\t%s\n" % (label, db))
        pct = extract_circular_rna_features(db)["pct_bulge"]
        score = (pct - 12.0) / 5.21
        exp.write("%s\t%d\t%.1f\t%.4f\n" % (label, len(db), pct, score))
PY

"$bin" --threads 1 validate_vleness/quick.tsv validate_vleness/quick_vleness.tsv
"$bin" --threads 4 validate_vleness/quick.tsv validate_vleness/quick_threads.tsv
cmp -s validate_vleness/quick_vleness.tsv validate_vleness/quick_threads.tsv
cmp -s validate_vleness/quick_vleness.tsv validate_vleness/quick_expected.tsv

# RNAfold: full label is cut at the first space, and the energy is stripped.
{
    printf '>bulge extra\n'
    printf 'ACGUACGU\n'
    printf '(.(...)) (-1.20)\n'
} > validate_vleness/rnafold.in
"$bin" --rnafold validate_vleness/rnafold.in validate_vleness/rnafold.tsv
awk 'NR==1 || $1=="bulge"' validate_vleness/quick_expected.tsv > validate_vleness/rnafold_expected.tsv
cmp -s validate_vleness/rnafold.tsv validate_vleness/rnafold_expected.tsv

# LSC has no header. The same circle is scored twice.
bulge_line=$(awk -F '\t' '$1=="bulge" { print; exit }' validate_vleness/quick_expected.tsv)
bulge_db=$(awk -F '\t' '$1=="bulge" { print $2; exit }' validate_vleness/quick.tsv)
{
    printf 'bulge\t-1.20\t%s\t-2.40\t%s\n' "$bulge_db" "$bulge_db"
} > validate_vleness/lsc.in
"$bin" --lsc validate_vleness/lsc.in validate_vleness/lsc.tsv
"$bin" validate_vleness/lsc.in validate_vleness/lsc_auto.tsv
cmp -s validate_vleness/lsc.tsv validate_vleness/lsc_auto.tsv
python3 - "$bulge_line" <<'PY'
import sys
label, length, pct, score = sys.argv[1].split("\t")
open("validate_vleness/lsc_expected.tsv", "w", newline="\n").write(
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (label, length, pct, pct, score, score, score)
)
PY
cmp -s validate_vleness/lsc.tsv validate_vleness/lsc_expected.tsv

# A pipe must work: the reader never seeks.
cat validate_vleness/quick.tsv | "$bin" /dev/stdin validate_vleness/stdin.tsv
cmp -s validate_vleness/stdin.tsv validate_vleness/quick_expected.tsv

# A first row that is not a dot-bracket is a header. A later bad structure
# is a warning, and --maxerror 0 stops the run after writing accepted rows.
printf 'name\tfold\nok\t((...))\n' > validate_vleness/header.tsv
"$bin" validate_vleness/header.tsv validate_vleness/header.out
cmp -s validate_vleness/header.out <(printf 'label\tlength\tpct_bulge\tS\nok\t7\t0.0\t-2.3033\n')
printf 'good\t((.))\nbad\t((.)\n' > validate_vleness/warn.tsv
set +e
"$bin" --maxerror 0 --log validate_vleness/warn.log \
    validate_vleness/warn.tsv validate_vleness/warn.out \
    > /dev/null 2> validate_vleness/warn.err
warn_status=$?
set -e
if [[ "$warn_status" -eq 0 ]]; then
    echo "maxerror 0 exited 0" >&2
    exit 1
fi
cmp -s validate_vleness/warn.out <(printf 'label\tlength\tpct_bulge\tS\ngood\t5\t0.0\t-2.3033\n')
grep -q $'^warning\tbad\t' validate_vleness/warn.err
grep -q $'^fatal\t' validate_vleness/warn.err
grep -q $'^start\t' validate_vleness/warn.log
grep -q $'^warning\tbad\t' validate_vleness/warn.log
grep -q $'^fatal\t' validate_vleness/warn.log
grep -q $'^finish\t' validate_vleness/warn.log
if "$bin" --rnafold --tsv 1,2 a b > /dev/null 2>&1; then
    echo "--rnafold and --tsv were accepted together" >&2
    exit 1
fi

awk -F '\t' '
{
    sub(/\r$/, "")
}
NR == FNR {
    if (!($0 in order)) {
        n += 1
        order[$0] = n
    }
    next
}
($1 in order) && !seen[$1]++ {
    rows[order[$1]] = $0
    found += 1
}
END {
    if (found != n) {
        printf "found %d of %d labels\n", found + 0, n + 0 > "/dev/stderr"
        exit 1
    }
    for (i = 1; i <= n; i++) {
        print rows[i]
    }
}
' big_output/control_labels_100_per_family.txt \
  big_output/controls_label_seq_ss.tsv \
  > validate_vleness/input.tsv

"$bin" --tsv 1,3 --threads 4 \
    validate_vleness/input.tsv \
    validate_vleness/vleness.tsv

python3 - <<'PY'
import sys
sys.path.insert(0, "py")
from dotbracket_circular_features import extract_circular_rna_features

with open("validate_vleness/input.tsv", "r", encoding="utf-8", newline="") as inp, \
     open("validate_vleness/python.tsv", "w", encoding="utf-8", newline="\n") as out:
    out.write("label\tlength\tpct_bulge\tS\n")
    for line in inp:
        line = line.rstrip("\n").rstrip("\r")
        if not line:
            continue
        label, _sequence, dot_bracket = line.split("\t")
        pct = extract_circular_rna_features(dot_bracket)["pct_bulge"]
        score = (pct - 12.0) / 5.21
        out.write("%s\t%d\t%.1f\t%.4f\n" % (label, len(dot_bracket), pct, score))
PY

if cmp -s validate_vleness/vleness.tsv validate_vleness/python.tsv; then
    rows=$(wc -l < validate_vleness/vleness.tsv)
    {
        echo "identical"
        echo "rows ${rows}"
    } > validate_vleness/summary.txt
else
    {
        echo "mismatch"
        diff -u validate_vleness/python.tsv validate_vleness/vleness.tsv | head -n 40
    } > validate_vleness/summary.txt
    exit 1
fi
