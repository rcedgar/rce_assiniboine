#!/usr/bin/env python3

# S       ROD|SRR25607789_680210_circle_680210:orf22:code1:1952-6820/5414(-1)     ROD|SRR25607789_680210_circle_680210:orf22:code1:1952-6820/5414(-1)     100
# H       ROD|SRR25607790_65861_circle_65861:orf18:code1:5104-9972/5416(-1)       ROD|SRR25607789_680210_circle_680210:orf22:code1:1952-6820/5414(-1)     100.00

centroid2members = {}

fn = "../big_output/orfs.id70.merged.tsv"
for line in open(fn):
	flds = line[:-1].split('\t')
	centroid = flds[2]
	member = flds[1]
	members = centroid2members.get(centroid, [])
	members.append(member)
	centroid2members[centroid] = members

cluster_nr = 0
for centroid, members in sorted(centroid2members.items(), key=lambda item: len(item[1]), reverse=True):
	cluster_nr += 1
	nmem = len(members)
	if nmem == 1:
		break
	group2n = {}
	nCe = 0
	nCem = 0
	nCo = 0
	nCom = 0
	nROD = 0
	nControl = 0
	for member in members:
		group = member.split('|')[0]
		if group == "ROD":
			nROD += 1
		else:
			nControl += 1
		n = group2n.get(group, 0)
		group2n[group] = n+1
		if member.endswith(")Ce"):
			nCe += 1
		elif member.endswith(")Cem"):
			nCem += 1
		elif member.endswith(")Co"):
			nCo += 1
		elif member.endswith(")Com"):
			nCom += 1
	s = "cluster_nr=%d" % cluster_nr
	s += "\tmembers=%d" % nmem
	s += "\tcentroid=%s" % centroid
	s += "\tgroups="
	for group, n in sorted(group2n.items(), key=lambda x: x[0], reverse=True):
		s += "%s:%d;" % (group, n)
	if nCe > 0:
		s += "\tnCe=%d(%.1f%%)" % (nCe, 100*nCe/nmem)
	if nCem > 0:
		s += "\tnCem=%d(%.1f%%)" % (nCem, 100*nCem/nmem)
	if nCo > 0:
		s += "\tnCo=%d(%.1f%%)" % (nCo, 100*nCo/nmem)
	if nCom > 0:
		s += "\tnCom=%d(%.1f%%)" % (nCom, 100*nCom/nmem)
	if nROD > 0 and nControl > 0:
		s += "\tnew_controls=%d" % nROD
	print(s)
