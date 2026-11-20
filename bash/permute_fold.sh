#!/usr/bin/env bash
# permute_fold.sh
#
# Shuffle a sequence N times (mononucleotide shuffle, composition preserved),
# fold the original and every permutation with `RNAfold --circ`, then print an
# ASCII histogram of the MFEs and the rank of the original among the permutations.
#
# Usage:  ./permute_fold.sh <input.fa> [n_permutations=1000] [seed=random]
#
# Outputs (next to the input file):
#   <name>.permuted.<N>.fa          permuted sequences
#   <name>.permuted.<N>.fold        raw RNAfold output (original first, then permutations)
#   <name>.permuted.<N>.mfe.tsv     id <tab> MFE (kcal/mol)
#   <name>.permuted.<N>.report.txt  the report that is also printed to stdout

set -euo pipefail

usage() { sed -n '2,14p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

[[ $# -ge 1 && $# -le 3 ]] || usage
IN=$1
N=${2:-1000}
SEED=${3:-$(( (RANDOM << 15) | RANDOM ))}
THREADS=${THREADS:-$(nproc 2>/dev/null || echo 1)}
NBINS=${NBINS:-20}      # approximate number of histogram bins
BARWIDTH=${BARWIDTH:-60} # width of the longest histogram bar

[[ -r $IN ]] || { echo "ERROR: cannot read '$IN'" >&2; exit 1; }
[[ $N =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: n_permutations must be a positive integer" >&2; exit 1; }
[[ $SEED =~ ^[0-9]+$ ]] || { echo "ERROR: seed must be a non-negative integer" >&2; exit 1; }
command -v RNAfold >/dev/null || { echo "ERROR: RNAfold not found in PATH" >&2; exit 1; }

DIR=$(dirname "$IN")
NAME=$(basename "$IN")
NAME=${NAME%.*}
PERM_FA="$DIR/$NAME.permuted.$N.fa"
FOLD_OUT="$DIR/$NAME.permuted.$N.fold"
MFE_TSV="$DIR/$NAME.permuted.$N.mfe.tsv"
REPORT="$DIR/$NAME.permuted.$N.report.txt"

# ---- read the input (first record only; multi-line sequence is joined) ------
NREC=$(grep -c '^>' "$IN" || true)
[[ $NREC -ge 1 ]] || { echo "ERROR: no FASTA header found in '$IN'" >&2; exit 1; }
[[ $NREC -eq 1 ]] || echo "WARNING: $NREC records in '$IN'; using the first only" >&2

SEQ_ID=$(awk '/^>/ { sub(/^>/, ""); print $1; exit }' "$IN")
SEQ=$(awk '/^>/ { n++; next } n == 1 { gsub(/[ \t\r]/, ""); printf "%s", $0 }' "$IN")
[[ -n $SEQ ]] || { echo "ERROR: empty sequence in '$IN'" >&2; exit 1; }

# ---- 1. permutations (seeded Fisher-Yates shuffle) --------------------------
awk -v seq="$SEQ" -v n="$N" -v seed="$SEED" -v id="$SEQ_ID" '
BEGIN {
    srand(seed)
    L = length(seq)
    w = length(n "")
    for (p = 1; p <= n; p++) {
        for (i = 1; i <= L; i++) a[i] = substr(seq, i, 1)
        for (i = L; i > 1; i--) {
            j = int(rand() * i) + 1
            t = a[i]; a[i] = a[j]; a[j] = t
        }
        printf ">%s_perm_%0*d\n", id, w, p
        for (i = 1; i <= L; i++) printf "%s", a[i]
        printf "\n"
    }
}' > "$PERM_FA"

# ---- 2. fold original + permutations as circular RNA ------------------------
{ printf '>%s\n%s\n' "$SEQ_ID" "$SEQ"; cat "$PERM_FA"; } \
    | RNAfold --circ --noPS --jobs="$THREADS" > "$FOLD_OUT"

# id <tab> MFE ; the structure line ends in "(-12.34)" or "( -2.34)"
awk '
/^>/ { sub(/^>/, ""); id = $1; next }
match($0, /\( *-?[0-9]+\.[0-9]+\)$/) {
    e = substr($0, RSTART, RLENGTH); gsub(/[() ]/, "", e)
    print id "\t" e
}' "$FOLD_OUT" > "$MFE_TSV"

NFOLDED=$(wc -l < "$MFE_TSV")
[[ $NFOLDED -eq $((N + 1)) ]] || {
    echo "ERROR: expected $((N + 1)) MFE values, parsed $NFOLDED from $FOLD_OUT" >&2; exit 1; }

ORIG_STRUCT=$(awk 'NR == 3 { print $1; exit }' "$FOLD_OUT")

# ---- 3. report: histogram + rank --------------------------------------------
awk -F'\t' -v nbins="$NBINS" -v barw="$BARWIDTH" -v id="$SEQ_ID" -v seed="$SEED" \
    -v infile="$IN" -v permfa="$PERM_FA" -v seq="$SEQ" -v struct="$ORIG_STRUCT" '
function floor(x,   i) { i = int(x); return (i > x) ? i - 1 : i }
function rep(c, k,   s) { s = ""; while (k-- > 0) s = s c; return s }
function fold60(s,   i) { for (i = 1; i <= length(s); i += 60) print "    " substr(s, i, 60) }
NR == 1 { orig = $2 + 0; next }
{ n++; e[n] = $2 + 0; sum += e[n]; sumsq += e[n] * e[n] }
END {
    # --- summary statistics of the permuted set
    mean = sum / n
    sd = (n > 1) ? sqrt((sumsq - n * mean * mean) / (n - 1)) : 0
    min = max = e[1]
    for (i = 1; i <= n; i++) {
        if (e[i] < min) min = e[i]
        if (e[i] > max) max = e[i]
        if (e[i] < orig) lower++
        else if (e[i] == orig) tied++
    }
    rank = lower + 1                      # rank 1 = lowest (most stable) MFE
    pval = (lower + tied + 1) / (n + 1)   # empirical P(permuted MFE <= original)

    # --- bins: "nice" width covering permuted values and the original
    lo = (orig < min) ? orig : min
    hi = (orig > max) ? orig : max
    raw = (hi - lo) / nbins
    if (raw <= 0) raw = 0.1
    mag = 10 ^ floor(log(raw) / log(10))
    f = raw / mag
    w = ((f <= 1) ? 1 : (f <= 2) ? 2 : (f <= 5) ? 5 : 10) * mag
    if (w < 0.1) w = 0.1
    start = floor(lo / w + 1e-9) * w
    nb = floor((hi - start) / w + 1e-9) + 1
    for (i = 1; i <= n; i++) cnt[floor((e[i] - start) / w + 1e-9)]++
    ob = floor((orig - start) / w + 1e-9)
    for (b = 0; b < nb; b++) if (cnt[b] > cmax) cmax = cnt[b]
    scale = (cmax > barw) ? cmax / barw : 1

    line = rep("=", 78)
    print line
    print "RNAfold --circ MFE: original vs. " n " permuted sequences"
    print line
    printf "Input file     : %s\n", infile
    printf "Sequence       : %s (%d nt)\n", id, length(seq)
    printf "Permutations   : %d, mononucleotide shuffle, seed %s\n", n, seed
    printf "Permuted fasta : %s\n\n", permfa

    print "Original sequence and circular MFE structure:"
    fold60(seq); print ""; fold60(struct); print ""

    printf "Original MFE   : %8.2f kcal/mol\n", orig
    printf "Permuted MFE   : mean %.2f, sd %.2f, min %.2f, max %.2f kcal/mol\n\n", mean, sd, min, max

    printf "Histogram of permuted MFEs (bin width %g kcal/mol, each # ~ %.2g sequences)\n\n", w, scale
    printf "  %19s  %6s\n", "MFE bin (kcal/mol)", "count"
    for (b = 0; b < nb; b++) {
        len = int(cnt[b] / scale + 0.5)
        if (cnt[b] > 0 && len == 0) len = 1
        printf "  [%8.2f,%8.2f)  %6d  %s%s\n", start + b * w, start + (b + 1) * w, cnt[b], rep("#", len), \
            (b == ob) ? sprintf("  <== ORIGINAL (%.2f)", orig) : ""
    }

    print ""
    print "Rank order (1 = lowest / most stable MFE, among original + permutations):"
    printf "  Original rank          : %d of %d\n", rank, n + 1
    printf "  Permuted with MFE <  original : %d (%.1f%%)\n", lower, 100 * lower / n
    printf "  Permuted with MFE == original : %d\n", tied
    printf "  Permuted with MFE >  original : %d (%.1f%%)\n", n - lower - tied, 100 * (n - lower - tied) / n
    if (sd > 0) printf "  Z-score                : %.2f\n", (orig - mean) / sd
    printf "  Empirical p (MFE <= original) : %.4g\n", pval
    print line
}' "$MFE_TSV" | tee "$REPORT"
