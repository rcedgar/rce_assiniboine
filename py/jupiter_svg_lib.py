#!/usr/bin/python3

import sys
import math

LINE_WIDTH = 1
BEND = 1
HANDLE_B = 1
COLORS = [ 0x0000ff, 0x00ff00, 0xffff00, 0xff0000 ]

# Dimensions are for canonical genome length, figure will be bigger
# or smaller in proportion to actual length
CIRCLE_DIAMETER = 300
CANONICAL_GENOME_LENGTH = 600
HORIZ_MARGIN = 5
VERT_MARGIN = 5
SCALE = 3
FLARE = 3
HeatMapValues =	[ 0, 0.333,	0.666, 1]

Title = None
if len(sys.argv) > 3:
	Title = sys.argv[3]

def ReadConfigFile(FileName):
	global LINE_WIDTH
	global BEND
	global HANDLE_B
	global COLORS

	for Line in open(FileName):
		Fields = Line[:-1].split("\t")
		Param = Fields[0]
		if Param == "line":
			LINE_WIDTH = float(Fields[1])
		elif Param == "bend":
			BEND = float(Fields[1])
		elif Param == "hb":
			HANDLE_B = float(Fields[1])
		elif Param == "flare":
			FLARE = float(Fields[1])
		elif Param == "colors":
			for i in range(0, 4):
				x = int(Fields[1+i], 0)
				COLORS[i] = x

if len(sys.argv) > 2:
	ConfigFileName = sys.argv[2]
	ReadConfigFile(ConfigFileName)

def HexToRGB(H):
	B = H%0x100
	h = (H - B)//0x100
	G = h%0x100
	h = h - G
	R = h//0x100

	assert R <= 255
	assert G <= 255
	assert B <= 255

	return R, G, B

Rs = []
Gs = []
Bs = []
for i in range(0, 4):
	R, G, B = HexToRGB(COLORS[i])
	Rs.append(R)
	Gs.append(G)
	Bs.append(B)

def GetSvgCircle(x, y, r, LineColor = "black", FillColor = "none"):
	return r'<circle cx="%.6g" cy="%.6g" r="%.6g" stroke-width="%.6g" stroke="%s" fill="%s" />' % \
	  (x, y, r, LINE_WIDTH, LineColor, FillColor)

def GetSvgLine(x1, y1, x2, y2, Color = "black"):
	return r'<line x1="%.6g" y1="%.6g" x2="%.6g" y2="%.6g" stroke-width="%.6g" stroke="%s" />' % \
	  (x1, y1, x2, y2, LINE_WIDTH, Color)

# https://medium.com/@bragg/cubic-bezier-curves-with-svg-paths-a326bb09616f
def GetBezierHandleCoords(Cx, Cy, Xi, Yi, Xj, Yj):
	Mx = (Xi + Xj)/2
	My = (Yi + Yj)/2

	Diff_x = Xi - Xj
	Diff_y = Yi - Yj
	Dist_ij = math.sqrt(Diff_x**2 + Diff_y**2)

	B0 = BEND*0.002

	dx = (Cx - Mx)*B0*Dist_ij
	dy = (Cy - My)*B0*Dist_ij

	Dx = Xi + dx
	Dy = Yi + dy

	Ex = Xj + dx
	Ey = Yj + dy

	MidDEx = (Dx + Ex)/2
	MidDEy = (Dy + Ey)/2

	B1 = HANDLE_B
	B2 = 1
	B12 = B1 + B2

	Fx = (B1*Dx + B2*MidDEx)/B12
	Fy = (B1*Dy + B2*MidDEy)/B12

	Gx = (B1*Ex + B2*MidDEx)/B12
	Gy = (B1*Ey + B2*MidDEy)/B12

	return Fx, Fy, Gx, Gy

######################################################
# GetColor(): Generate heat map color by interpolation
######################################################
def GetColor(Value):
	assert Value >= 0 and Value <= 1
	N = len(HeatMapValues)
	assert len(Rs) == N
	assert len(Gs) == N
	assert len(Bs) == N

	for i in range(1, N):
		LoValue = HeatMapValues[i-1]
		HiValue = HeatMapValues[i]
		assert LoValue < HiValue
		if Value >= LoValue and Value <= HiValue:
			f = float(Value - LoValue)/float(HiValue - LoValue)
			assert f >= 0.0 and f <= 1.0
			R = int((Rs[i] - Rs[i-1])*f + Rs[i-1])
			G = int((Gs[i] - Gs[i-1])*f + Gs[i-1])
			B = int((Bs[i] - Bs[i-1])*f + Bs[i-1])
			assert R >= 0 and R <= 255
			assert G >= 0 and G <= 255
			assert B >= 0 and B <= 255
			Str = "#%02x%02x%02x" % (R, G, B)
			return Str
	assert False
######################################################

def fix_rounding(x):
	if x < 0:
		return 0
	if x > 1:
		return 1
	return x

def ReadTSVFile(TSVFileName):
	global Xs, Ys, Ss, Map_i, j, Pairs, IsLoops

	f = open(TSVFileName)
	Header = f.readline()
	assert Header.startswith("i\tx\ty\tS\tj")

	Xs = []
	Ys = []
	Ss = []
	Map_i_j = {}
	Pairs = []
	IsLoops = []

	Counter = 0
	while 1:
		Line = f.readline()
		if len(Line) == 0:
			break
		Fields = Line[:-1].split('\t')
		assert len(Fields) >= 5
	
		i = int(Fields[0])
		x = float(Fields[1])
		y = float(Fields[2])
		S = float(Fields[3])
		sj = Fields[4]
		IsLoop = None
		if sj == ".":
			j = None
		else:
			j = int(sj)
			Map_i_j[i] = j
			Map_i_j[j] = i
	
			if i < j:
				Pair = (i, j)
			else:
				Pair = (j, i)
			Pairs.append(Pair)

			if len(Fields) > 5:
				strIsLoop = Fields[5]
				if strIsLoop == "Rod":
					IsLoop = False
				elif strIsLoop == "Loop":
					IsLoop = True
				else:
					assert False

		assert i == Counter
		Counter += 1
		x = fix_rounding(x)
		y = fix_rounding(y)

		Xs.append(x)
		Ys.append(y)
		Ss.append(S)
		IsLoops.append(IsLoop)

def DrawJupiter(TSVFileName, JupiterCenter_x, JupiterCenter_y, JupiterRadius):

	print((r'<circle cx="%.6g" cy="%.6g" r="%.6g" stroke-width="%.6g" stroke="%s" fill="%s" />' % \
	  (JupiterCenter_x, JupiterCenter_y, JupiterRadius, 3.75, "black", "white")))
	ReadTSVFile(TSVFileName)

	print((r'<g id="%s">' % TSVFileName))

	PairCount = len(Pairs)
	for PairIndex in range(0, PairCount):
		Pair = Pairs[PairIndex]
		i, j = Pair
		Xi = Xs[i]
		Yi = Ys[i]
		Xj = Xs[j]
		Yj = Ys[j]
		Si = Ss[i]
		Sj = Ss[j]
		IsLoop_i = IsLoops[i]
		IsLoop_j = IsLoops[j]
		IsLoop = IsLoop_i is True or IsLoop_j is True

		assert Si == Sj

		Diam = 2*JupiterRadius
		PLOT_Xi = JupiterCenter_x + (Xi - 0.5)*Diam
		PLOT_Yi = JupiterCenter_y + (Yi - 0.5)*Diam

		PLOT_Xj = JupiterCenter_x + (Xj - 0.5)*Diam
		PLOT_Yj = JupiterCenter_y + (Yj - 0.5)*Diam

		if Si < 0:
			Si = 0
		if Si > 1:
			Si = 1
		Color = GetColor(Si)
	
		Handle_ix, Handle_iy, Handle_jx, Handle_jy = \
		  GetBezierHandleCoords(JupiterCenter_x, JupiterCenter_y, PLOT_Xi, PLOT_Yi, PLOT_Xj, PLOT_Yj)

		print((r'<path d="M%.6g,%.6g C%.6g,%.6g %.6g,%.6g %.6g,%.6g" stroke-width="%.6g" fill="none" stroke="%s" />' \
			% (PLOT_Xi, PLOT_Yi, Handle_ix, Handle_iy, Handle_jx, Handle_jy, PLOT_Xj, PLOT_Yj, LINE_WIDTH, Color)))

	#	print(r'<line x1="%.6g" y1="%.6g" x2="%.6g" y2="%.6g" stroke-width="%.6g" stroke="%s" />' % \
	#	  (PLOT_Xi, PLOT_Yi, PLOT_Xj, PLOT_Yj, 1, "red"))
	print("</g>")
