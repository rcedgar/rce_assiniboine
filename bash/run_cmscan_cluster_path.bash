#!/bin/bash -e

path=$1

cm=~/rce_assiniboine/data/small-2026-11-06.cm

echo $path

fa=$path/circles.maxlen5000.rotated.pctid90.fa

cmscan \
  -o /dev/null \
  -E 10 \
  --tblout $path/cmscan.tbl \
  $cm $fa

python ~/rce_assiniboine/py/cmscan_tbl2fev.py \
	$path/cmscan.tbl 10 \
	> $path/cmscan.fev

ls -lh $path/cmscan.*
