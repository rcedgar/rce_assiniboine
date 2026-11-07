#!/usr/bin/env python3

# Construct directories e.g. 
#   rce_assiniboine/big_clusters/2026-11-06/nn
# One directory per cluster in fevfn
# Directory contains cirlces.fa and orfs.fa

import sys
import fasta
from pathlib import Path

# fev
# cluster_nr=1128 members=91      centroid=plasmidsub7000_gb260921|NC_008352.1:orf1:code1:91-1092/2278(+1)  ...

# tsv
# S       ROD|SRR23617204_519500_circle_519500:orf24:code1:656-5803/5539(-1)      ROD|SRR23617204_519500_circle_519500:orf24:code1:656-5803/5539(-1)      100
# H       ROD|SRR22198339_1843530_circle_1843530:orf22:code1:1057-6204/5539(-1)   ROD|SRR23617204_519500_circle_519500:orf24:code1:656-5803/5539(-1)      100.00

# nt
# >ROD|DRR000346_804025_circle_804025
# GGTCGTGCCACGCAGGTAATCACCGCCAGCCCGTGATCCGATCACCTTTTGAACACCGGAGACTGTTTCGCTGCCATAGGCACCGGTGACGGTCCCGGCCGACAAGTTGACTTGTACAGCACCGGTCGCATTACCGAAGGTGAGCGTGTTATTCGTGCCCGCACCAGTAATGGTAATAGCAGTGGCCCCACCGCCATCGAAGGTCACGTTGCCAGCGGATCCCTGTAACAAGGTCGTATT
# >ROD|DRR000426_1011226_circle_1011226
# GCACGTGTGCATACAGCATTGTGTGTGTACATGCGTGTGCGTGCGCATACAGTATTGTGTGTGTGTACACGTGCGCACGTGCGCATACAGCATTGTGTGTGTACATGTGCGTGCACGTGCTCATACAGCATTGTGTGCGTACATGCGCGCGCGTGTGCATACAGCATTGTGTGTGTGTACCTGCGCGCATGTGCGCATACAGCATTGTGTGTGTACATGCACGCGCGTGCGCATACAGTATTGTGTGTGTGTGTACATGCGT

# orf
# >ROD|DRR235453_664787_circle_664787:orf1:code1:3-329/327(+3)Cem
# NNVSTPANEVITVFNVQEGPTATPSSSSGNEDSDIGVSLSGTDPDGTVVSVTVTSLPPASQGVLYLADGTTPVVAGVALT
# AAQAASLVFSPAADFHGSVTIPFTVTDDE

fevfn = sys.argv[1]     # rce_assiniboine/analysis/2026-11-06/top200_selected.fev
tsvfn = sys.argv[2]     # rce_assiniboine/big_output/big_output/orfs.id70.merged.tsv
ntfn = sys.argv[3]      # rce_assiniboine/big_output/all.fa
orffn = sys.argv[4]     # rce_assiniboine/big_output/dust_filtered.reorf_longest.fa
rootdir = sys.argv[5]   # rce_assiniboine/big_clusters/2026-11-06/

if not rootdir.endswith('/'):
	rootdir += '/'

centroids = set()
centroid2cluster_nr = {}

# Read fev
sys.stderr.write(fevfn + '\n')
for line in open(fevfn):
	flds = line[:-1].split('\t')
	assert flds[0].startswith("cluster_nr=")
	assert flds[2].startswith("centroid=")
	cluster_nr = int(flds[0].replace("cluster_nr=", ""))
	centroid = flds[2].replace("centroid=", "")
	centroid = centroid.split(":")[0]
	centroids.add(centroid)
	centroid2cluster_nr[centroid] = cluster_nr
	cluster_dir = rootdir + str(cluster_nr)
	Path(cluster_dir).mkdir(parents=True, exist_ok=True)

sys.stderr.write("%d selected clusters\n" % len(centroids))

members = set()
member2centroid = {}
sys.stderr.write(tsvfn + '\n')
nmem = 0
for line in open(tsvfn):
	flds = line[:-1].split('\t')
	assert len(flds) == 4
	assert flds[0] == 'S' or flds[0] == 'H'
	centroid = flds[2]
	centroid = centroid.split(":")[0]
	if centroid not in centroids:
		continue
	nmem += 1
	member = flds[1]
	member = member.split(":")[0]
	members.add(member)
	member2centroid[member] = centroid
sys.stderr.write("%d members\n" % nmem)

n_nt_seqs = 0
def onseq_nt(label, seq):
	global n_nt_seqs
	if not label in members:
		return
	n_nt_seqs += 1
	member = label
	centroid = member2centroid[member]
	cluster_nr = centroid2cluster_nr[centroid]
	cluster_dir = rootdir + str(cluster_nr)
	nt_fasta_fn = cluster_dir + "/circles.fa"
	f = open(nt_fasta_fn, "a")
	fasta.WriteSeq(f, seq, label)
	f.close()

sys.stderr.write(ntfn + '\n')
fasta.ReadSeqsOnSeq(ntfn, onseq_nt)
sys.stderr.write("%d circles\n" % n_nt_seqs)

n_orf_seqs = 0
def onseq_orf(label, seq):
	global n_orf_seqs
	member = label.split(":")[0]
	if not member in members:
		return
	n_orf_seqs += 1
	centroid = member2centroid[member]
	cluster_nr = centroid2cluster_nr[centroid]
	cluster_dir = rootdir + str(cluster_nr)
	orf_fasta_fn = cluster_dir + "/orfs.fa"
	f = open(orf_fasta_fn, "a")
	fasta.WriteSeq(f, seq, label)
	f.close()

sys.stderr.write(orffn + '\n')
fasta.ReadSeqsOnSeq(orffn, onseq_orf)
sys.stderr.write("%d orf seqs\n" % n_orf_seqs)

