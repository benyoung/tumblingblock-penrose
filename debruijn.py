#!/usr/bin/python
# Construct a tiling which looks like tumbling block on the left hand side, and
# like a penrose tiling on the right hand side.
#
# Usage:  debruijn.py > tumblingblockpenrose.eps
#
# This also creates (as a side effect) debruijn.eps, a picture of the de Bruijn
# lines from which the tiling is derived.

from subprocess import call
import sys
import math
import cmath
import copy
from collections import defaultdict

# Set up slopes and y-intercepts of De Bruijn lines.
dim = 5
theta = cmath.pi / dim
I = complex(0,1)
xi = cmath.exp(I * theta)
def slope(z):
    return z.imag / z.real

slopes = [slope(xi**k) for k in range(dim)]
colors = ["0 0 0", "1 0 0", "0 0 0", "0 0 0", "0 0 1"]

def webcolor(color):
    hexcolors = [str(color[0:2]), str(color[2:4]), str(color[4:6])]
    (r,g,b) = [(int(x,16) + 0.0) / 255.0 for x in hexcolors]
    return "%0.3f %0.3f %0.3f"%(r,g,b)
    
tilecolors = {
# tumbling block tile
    frozenset((0,2)):webcolor("0D0D01"), 
    frozenset((0,3)):webcolor("222214"), 
    frozenset((2,3)):webcolor("47472C"),

# xi lines
    frozenset((0,4)): webcolor("856AD9"),
    frozenset((3,4)): webcolor("1B0072"),
    frozenset((2,4)): webcolor("2900AD"),

# xi^4 lines
    frozenset((0,1)): webcolor("FD7F7C"),
    frozenset((1,2)): webcolor("A60400"),
    frozenset((1,3)): webcolor("FB0600"),

# (xi, xi^4) tiles
    frozenset((1,4)): webcolor("00E900"),
}


debruijn = []

# Penrose tiling.
increment = [1/((xi**k).real) for k in range(dim)]
normal_intercepts = [(I * xi**k).imag for k in range(dim) ]
normal_offsets = [0.5, -1.55, 0.3, 0.4, 0.35, ] # should sum to 0
offset = [normal_offsets[k] * normal_intercepts[k] for k in range(dim)]
indices = range(-8,8)

debruijn.append([offset[0] + u for u in indices])                           # direction 1
debruijn.append([offset[1] + increment[1]*t for t in range(-4,6)])    # direction xi
debruijn.append([offset[2] + increment[2]*t for t in range(-6,18)])    # direction xi^2
debruijn.append([offset[3] + increment[3]*t for t in range(-5,18)])    # direction xi^3
debruijn.append([offset[4] + increment[4]*t for t in range(-4,5)])    # direction xi^4

# Let's add some tumbling block stuff on there.
tb_offset2 = debruijn[2][-1]
tb_offset3 = debruijn[3][-1]

debruijn[2].extend([tb_offset2 + 2*t for t in range(0,15)])
debruijn[3].extend([tb_offset3 - 2*t for t in range(0,15)])


for k in range(1,dim):
    intercepts = [k*0.1 + 1/((xi**k).real) * t for t in indices]
    debruijn.append(intercepts) 



# X coordinate of intersection of two lines.  Each line is an ordered pair (m,b)
def intersect(L1, L2):
    return (L2[1]-L1[1])/(L1[0]-L2[0])

def inbox(p, bb):
    return ((bb[0] < p[0]) and (bb[1] < p[1]) and (bb[2] > p[0]) and (bb[3] > p[1]))

# xmin, ymin, xmax, ymax
eps = 0.00000001
bb = (-300, offset[0] + indices[0] -eps, 300, offset[0] + indices[-1]  +eps)

# Render De Bruijn lines.
def debruijn_render(debruijn, bb, f):
    (xmin, ymin, xmax, ymax) = bb
    print >> f, "%!PS-Adobe-2.0 EPSF-2.0"
    print >> f, "%%BoundingBox: 122 368 496 435" 
    print >> f, "500 400 translate"
    print >> f, "15 15 scale"
    print >> f, "0.01 setlinewidth"

    for k in range(dim): 
        print >> f, "%s setrgbcolor" % colors[k]
        m = slopes[k]
        for b in debruijn[k]:
            x0 = xmin + eps
            y0 = m*x0 + b
            if not inbox((x0, y0), bb):
                y0 = ymin + eps
                x0 = (y0 - b) / m
                if not inbox((x0, y0), bb):
                    y0 = ymax - eps
                    x0 = (y0 - b) / m

            x1 = xmax - eps
            y1 = m*x1 + b
            if not inbox((x1, y1), bb):
                y1 = ymax - eps
                x1 = (y1 - b) / m
                if not inbox((x1, y1), bb):
                    y1 = ymin + eps
                    x1 = (y1 - b) / m
            print >> f, "newpath %0.3f %0.3f moveto %0.3f %0.3f lineto closepath stroke" %(x0, y0, x1, y1)

# Find all pairwise intersections between De Bruijn lines.
# return a dictionary of pairs:
#
# debruijn line -> ( intersection point, intersecting debruijn line) 
#
# where a debruijn line is specified by an ordered pair 
#     (index into slopes array, index into y-intercepts array)
def find_intersections(debruijn, slopes, bb):
    intersections = defaultdict(list)

    for k1 in range(dim):
        m1 = slopes[k1]
        for j1 in range(len(debruijn[k1])):
            b1 = debruijn[k1][j1]
            for k2 in range(k1+1, dim):
                m2 = slopes[k2]
                for j2 in range(len(debruijn[k2])):
                    b2 = debruijn[k2][j2]
                    x = intersect((m1,b1), (m2,b2))
                    y = m1 * x + b1
                    if inbox((x,y), bb):
                        intersections[(k1, j1)].append(((x, y), (k2, j2)))
                        intersections[(k2, j2)].append(((x, y), (k1, j1)))
    return intersections

def show_intersections(intersections, f):
    radius = 0.03
    for (k,j) in intersections:
        line = intersections[(k,j)]
        for record in line:
            #radius += 0.005
            p=record[0]
            k2=record[1][0]
            color = tilecolors[frozenset((k,k2))]

            print >> f, "%s setrgbcolor" % color
            print >> f, "newpath %0.3f %0.3f %0.3f 0 360 arc fill" % (p[0], p[1], radius)

def sort_intersections(intersections):
    for v in intersections.values():
        v.sort()

# We only needed to compute (x,y) coordinates in the debruijn lines to render them,
# and to rank them.  Once that's done, we don't need them anymore.  
# 
# When this function is done, intersections is just a dictionary of lists of lines L.
def discard_debruijn_coordinates(intersections):
    for L in intersections:
        intersections[L] = [x[1] for x in intersections[L]]


# Place tiles cooresponding to intersections on a de bruijn line.
# Starting point is:
#    lower left corner of leftmost tile for nonnegative slope.
#    lower right corner of rightmost tile for negative slope.
def place_tiles_along_line(start, L, slopes, crossings, coords):

    m = slopes[L[0]]
    v = I * complex(1, m) # positive normal to (direction vector for L with positive x coord).   
    v = v / abs(v)        # make that a unit normal

    direction = 1
    if(m < -eps):  # if slope is negative, do everything backwards.
        direction = -1
        crossings.reverse()

    e1 = (start, start + v)

    for L2 in crossings:
        m2 = slopes[L2[0]]
        w = I * complex(1,m2) #again a unit normal to L2
        w = w / abs(w)
        # Is w the right orientation?  Suppose first direction = 1.  
        # Then w should be on the "positive" side of v.
        # So we should have Arg v - pi < Arg w < Arg v.
        # That's the same as 0 < Arg v - Arg w, because of where the branch cut of arg is. 
        # 
        # If direction = -1, though, then w should be on the negative side of v.
        if 0 > direction * (cmath.phase(v) - cmath.phase(w)):
            w = -w

        e2 = (e1[0] + w, e1[1] + w)
        coords[(L, L2)] = (e1, e2)
        e1 = e2 

    


def render_tiles(intersection_coords):
    tiles = {}
    for I in intersection_coords:  # make intersections unordered
        slopes = frozenset([L[0] for L in I])
        tiles[frozenset(I)] = (intersection_coords[I],  tilecolors[slopes])
    for (((p,q),(s,r)),color) in tiles.values():
        print >> f, "%s setrgbcolor" % color
        print >> f, "newpath %0.3f %0.3f moveto %0.3f %0.3f lineto"     % (p.real, p.imag, q.real, q.imag)
        print >> f, "  %0.3f %0.3f lineto %0.3f %0.3f lineto fill"%(r.real,r.imag,s.real,s.imag)
        print >> f, "0 setgray"
        print >> f, "newpath %0.3f %0.3f moveto %0.3f %0.3f lineto"     % (p.real, p.imag, q.real, q.imag)
        print >> f, "  %0.3f %0.3f lineto %0.3f %0.3f lineto closepath stroke"%(r.real,r.imag,s.real,s.imag)

intersections = find_intersections(debruijn, slopes, bb)
sort_intersections(intersections)
debruijnfile = open("debruijn.eps", 'w')
debruijn_render(debruijn, bb, debruijnfile)
show_intersections(intersections, debruijnfile)
print >> debruijnfile, "showpage"
debruijnfile.close()

discard_debruijn_coordinates(intersections)


# this header needs a LOT of tweaking to get a nice looking crop.
f = open("tumblingblockpenrose.eps", "w")
print >> f, "%!PS-Adobe-2.0 EPSF-2.0"
print >> f, "%%BoundingBox: 120 393 720 555" 

#debug
print >> f, "10 400 translate"
print >> f, "10 10 scale"
print >> f, "0.01 setlinewidth"

coords = {}

L = (0,0)
start = complex(0,0)
place_tiles_along_line(start, L, slopes, intersections[L], coords)

for L2 in intersections[L]:
    m = L2[0]
    tile = coords[(L, L2)]
    if(slopes[m] > 0):
        start = tile[1][0]
    else:
        start = tile[0][0]
    place_tiles_along_line(start, L2, slopes, intersections[L2], coords)

render_tiles(coords)

print >> f, "showpage"

f.close()

call("epspdf debruijn.eps", shell=True)
call("epspdf tumblingblockpenrose.eps", shell=True)



