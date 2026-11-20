#!/usr/bin/python3

import sys
import math
import os
import shutil
import tempfile
import argparse

# Input is tsv with optional header:
#	sample_name	length	mfe	dot_bracket
# Output is an SVG contact sheet of Jupiter plots.
# Every circle is drawn at the same --diameter. Length only changes the base spacing.
# argv is hidden before import: jupiter_svg_lib reads argv[2] as a config file.

parser = argparse.ArgumentParser(description="Jupiter plots of dot-bracket structures")
parser.add_argument("tsv", help="sample_name, length, mfe, dot_bracket")
parser.add_argument("--cols", type=int, default=4, help="plots per row")
parser.add_argument("--margin", type=float, default=0.05, help="margin on each side, as a fraction of diameter")
parser.add_argument("--diameter", type=float, default=None, help="circle diameter for every plot (default: CIRCLE_DIAMETER)")
Args = parser.parse_args()

InFile = Args.tsv
sys.argv = [ sys.argv[0] ]
import jupiter_svg_lib
import rnafold_dotbracket_rod

if Args.diameter == None:
	Args.diameter = jupiter_svg_lib.CIRCLE_DIAMETER
assert Args.cols >= 1
assert Args.margin >= 0
assert Args.diameter > 0

NCOL = Args.cols
STROKE = 3.75
FONT_SIZE = 12
LINE_H = 16
LABEL_H = 2*LINE_H
MARGIN_FRAC = Args.margin
TWOPI = 2.0*3.1415927

def XmlText(S):
	S = S.replace("&", "&amp;")
	S = S.replace("<", "&lt;")
	S = S.replace(">", "&gt;")
	return S

def WritePlotTSV(Path, Parens):
	L = len(Parens)
	Stack = []
	Partner = [ None ] * L
	Paired = [ False ] * L
	for i in range(0, L):
		c = Parens[i]
		if c == '(':
			Stack.append(i)
		elif c == ')':
			assert len(Stack) > 0
			i2 = Stack.pop()
			Partner[i2] = i
			Paired[i2] = True
			Paired[i] = True
		elif c == '.':
			continue
		else:
			assert False
	assert len(Stack) == 0

	# Both bases of a pair have S=1. DrawJupiter requires Ss[i] == Ss[j].
	# Only the opening base records j, so each pair is drawn once.
	f = open(Path, "w", newline="\n")
	f.write("i\tx\ty\tS\tj\n")
	for i in range(0, L):
		Theta = TWOPI*float(i)/float(L)
		x = (1.0 + math.cos(Theta))/2.0
		y = (1.0 + math.sin(Theta))/2.0
		j = Partner[i]
		if Paired[i]:
			S = 1
		else:
			S = 0
		if j == None:
			f.write("%d\t%.4f\t%.4f\t%d\t.\n" % (i, x, y, S))
		else:
			f.write("%d\t%.4f\t%.4f\t%d\t%d\n" % (i, x, y, S, j))
	f.close()

def EmitText(X, Y, SlotW, Text):
	# Keep a gap so neighboring titles do not meet.
	FitW = SlotW*0.92
	TextLen = ""
	if len(Text)*0.6*FONT_SIZE > FitW:
		TextLen = ' textLength="%.6g" lengthAdjust="spacingAndGlyphs"' % FitW
	print('<text x="%.6g" y="%.6g" text-anchor="middle" font-family="sans-serif" font-size="%d"%s>%s</text>' %
		(X, Y, FONT_SIZE, TextLen, Text))

Samples = []
SeenData = False
for Line in open(InFile):
	Line = Line.strip()
	if len(Line) == 0:
		continue
	Fields = Line.split('\t')
	if not SeenData:
		try:
			int(Fields[1])
		except ValueError:
			continue
	SeenData = True

	assert len(Fields) == 4
	Name = Fields[0]
	Length = int(Fields[1])
	Mfe = Fields[2]
	Parens = Fields[3]
	assert Length == len(Parens)
	Samples.append((Name, Length, Mfe, Parens))

assert len(Samples) > 0

# Every circle has the same diameter. Margin is that fraction of the diameter.
Radius = Args.diameter/2.0
Diam = Args.diameter
Margin = MARGIN_FRAC*Diam
SlotW = Diam + STROKE + 2.0*Margin
SlotH = LABEL_H + Diam + STROKE + 2.0*Margin
NRows = (len(Samples) + NCOL - 1)//NCOL
Width = NCOL*SlotW
Height = NRows*SlotH

# Highest rodness first. Equal scores keep input order.
Scored = []
for Sample in Samples:
	Name, Length, Mfe, Parens = Sample
	Score = rnafold_dotbracket_rod.Rodness(Parens)[0]
	Scored.append((Name, Length, Mfe, Parens, Score))
Scored.sort(key=lambda Item: -Item[4])

Plots = []
for Index in range(0, len(Scored)):
	Name, Length, Mfe, Parens, Score = Scored[Index]
	Col = Index % NCOL
	Row = Index // NCOL
	Cx = Col*SlotW + Margin + STROKE/2.0 + Radius
	Cy = Row*SlotH + LABEL_H + Margin + STROKE/2.0 + Radius
	LabelY1 = Row*SlotH + FONT_SIZE
	LabelY2 = Row*SlotH + FONT_SIZE + LINE_H
	Plots.append((Name, Length, Mfe, Parens, Score, Cx, Cy, Radius, LabelY1, LabelY2, SlotW))

print('<?xml version="1.0" encoding="UTF-8"?>')
print('<svg xmlns="http://www.w3.org/2000/svg" width="%.6f" height="%.6f" viewBox="0 0 %.6f %.6f">' %
	(Width, Height, Width, Height))
print('<rect width="100%" height="100%" fill="#ffffff"/>')

TmpDir = tempfile.mkdtemp(prefix="jupiter_")
try:
	for Index in range(0, len(Plots)):
		Name, Length, Mfe, Parens, Score, Cx, Cy, Radius, LabelY1, LabelY2, SlotW = Plots[Index]
		EmitText(Cx, LabelY1, SlotW, XmlText(Name))
		Stats = "L=%d MFE=%s Rodness=%.4f" % (Length, XmlText(Mfe), Score)
		EmitText(Cx, LabelY2, SlotW, Stats)
		PlotPath = os.path.join(TmpDir, "%d.tsv" % Index)
		WritePlotTSV(PlotPath, Parens)
		jupiter_svg_lib.DrawJupiter(PlotPath, Cx, Cy, Radius)
finally:
	shutil.rmtree(TmpDir)

print('</svg>')
