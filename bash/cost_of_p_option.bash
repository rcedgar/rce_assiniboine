#!/bin/bash -e

# Time RNAfold --circ --noPS --jobs=1 on the first 10 controls, with and without -p.

root="$(cd "$(dirname "$0")/.." && pwd)"
src="$root/big_output/controls.subset1000.fa"
outdir="$root/cost_of_p_option"
sample="$outdir/sample10.fa"
fallback_rnafold=/home/bob/miniconda3/envs/std/bin/RNAfold

if [ ! -s "$src" ]; then
	echo "ERROR: not found $src" >&2
	exit 1
fi

if ! command -v RNAfold >/dev/null 2>&1; then
	if [ ! -x "$fallback_rnafold" ]; then
		echo "ERROR: RNAfold not found in PATH or at $fallback_rnafold" >&2
		exit 1
	fi
	export PATH="$(dirname "$fallback_rnafold"):$PATH"
fi

mkdir -p "$outdir"

# First 10 FASTA records, in file order.
awk 'BEGIN { n = 0 } /^>/ { n++ } n <= 10 { print } n > 10 { exit }' "$src" > "$sample"

nseq=$(grep -c '^>' "$sample")
if [ "$nseq" -ne 10 ]; then
	echo "ERROR: expected 10 sequences, got $nseq" >&2
	exit 1
fi

# Lengths of the 10 sampled sequences (joined multi-line records).
lengths=$(awk '
	/^>/ {
		if (started) print len
		started = 1
		len = 0
		next
	}
	{
		gsub(/[ \t\r]/, "")
		len += length($0)
	}
	END { if (started) print len }
' "$sample" | paste -sd, -)
min_len=$(echo "$lengths" | tr ',' '\n' | sort -n | head -n 1)
max_len=$(echo "$lengths" | tr ',' '\n' | sort -n | tail -n 1)

cd "$outdir"

/bin/time -v -o time_without_p.txt -- \
	RNAfold --circ --noPS --jobs=1 sample10.fa \
	> fold_without_p.txt 2> fold_without_p.err

/bin/time -v -o time_with_p.txt -- \
	RNAfold --circ --noPS --jobs=1 -p sample10.fa \
	> fold_with_p.txt 2> fold_with_p.err

# GNU time -v fields. Elapsed is h:mm:ss or m:ss; convert to seconds.
time_field() {
	local file=$1 label=$2
	awk -v label="$label" '
		{
			line = $0
			sub(/^[ \t]+/, "", line)
			prefix = label ": "
			if (index(line, prefix) == 1) {
				print substr(line, length(prefix) + 1)
				exit
			}
		}
	' "$file"
}

elapsed_seconds() {
	local raw=$1
	awk -v t="$raw" 'BEGIN {
		n = split(t, a, ":")
		if (n == 3) sec = a[1] * 3600 + a[2] * 60 + a[3]
		else if (n == 2) sec = a[1] * 60 + a[2]
		else sec = a[1]
		printf "%.2f", sec
	}'
}

wall_without=$(elapsed_seconds "$(time_field time_without_p.txt "Elapsed (wall clock) time (h:mm:ss or m:ss)")")
wall_with=$(elapsed_seconds "$(time_field time_with_p.txt "Elapsed (wall clock) time (h:mm:ss or m:ss)")")
user_without=$(time_field time_without_p.txt "User time (seconds)")
user_with=$(time_field time_with_p.txt "User time (seconds)")
sys_without=$(time_field time_without_p.txt "System time (seconds)")
sys_with=$(time_field time_with_p.txt "System time (seconds)")
rss_without=$(time_field time_without_p.txt "Maximum resident set size (kbytes)")
rss_with=$(time_field time_with_p.txt "Maximum resident set size (kbytes)")

awk -v nseq="$nseq" -v lengths="$lengths" -v min_len="$min_len" -v max_len="$max_len" \
	-v wall_without="$wall_without" -v wall_with="$wall_with" \
	-v user_without="$user_without" -v user_with="$user_with" \
	-v sys_without="$sys_without" -v sys_with="$sys_with" \
	-v rss_without="$rss_without" -v rss_with="$rss_with" \
	-v rnafold="$(command -v RNAfold)" \
'BEGIN {
	wall_delta = wall_with - wall_without
	rss_delta = rss_with - rss_without
	wall_ratio = (wall_without > 0) ? wall_with / wall_without : 0
	rss_ratio = (rss_without > 0) ? rss_with / rss_without : 0
	printf "RNAfold -p incremental cost\n"
	printf "RNAfold: %s\n", rnafold
	printf "input: first %d sequences of big_output/controls.subset1000.fa\n", nseq
	printf "lengths (nt): %s\n", lengths
	printf "length range (nt): %s-%s\n", min_len, max_len
	printf "command: RNAfold --circ --noPS --jobs=1 sample10.fa\n"
	printf "\n"
	printf "%-12s %12s %12s %12s %14s\n", "run", "wall_s", "user_s", "sys_s", "max_rss_kb"
	printf "%-12s %12s %12s %12s %14s\n", "without -p", wall_without, user_without, sys_without, rss_without
	printf "%-12s %12s %12s %12s %14s\n", "with -p", wall_with, user_with, sys_with, rss_with
	printf "\n"
	printf "wall time: -p adds %.2f s (%.2fx)\n", wall_delta, wall_ratio
	printf "max RSS:   -p adds %d kB (%.2fx)\n", rss_delta, rss_ratio
}' > report.txt

cat report.txt
