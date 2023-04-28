import pickle
import random
import sys

from shapely import affinity
from shapely.geometry import Point

# -----------------------------------------------------------------------------
inputfile = open(sys.argv[1], "r")

parameters = []
for parameter in inputfile.readlines():
    parameters.append(parameter.split(":")[1].strip(" \n"))

inputfile.close()

seed = int(parameters[0])
iterations = int(parameters[1])
rows = int(parameters[2])
columns = int(parameters[3])
tetrominolist = [tetro.strip(" ").upper() for tetro in parameters[4].split(",")]
anglelist = [int(text.strip(" ")) for text in parameters[5].split(",")]

idfile = parameters[6].split(",")[0]
idfields = [text.strip(" ").upper() for text in parameters[6].split(",", 1)[1].strip(" []").split(",")]
xyfile = parameters[7].split(",")[0]
xyfields = [text.strip(" ").upper() for text in parameters[7].split(",", 1)[1].strip(" []").split(",")]

pickledoutput = parameters[8]
countoutput = parameters[9]
foldsoutput = parameters[10]

# -----------------------------------------------------------------------------
random.seed(seed)

blank = [[0 for column in range(columns)] for row in range(rows)]

a = [(0, 0), (1, 0), (2, 0), (3, 0)]
b = [(0, 0), (0, 1), (0, 2), (0, 3)]

c = [(0, 0), (1, 0), (2, 0), (2, 1)]
d = [(0, 0), (1, 0), (1, -1), (1, -2)]
e = [(0, 0), (0, 1), (1, 1), (2, 1)]
f = [(0, 0), (0, 1), (0, 2), (1, 0)]

g = [(0, 0), (1, 0), (2, 0), (2, -1)]
h = [(0, 0), (1, 0), (1, 1), (1, 2)]
i = [(0, 0), (0, 1), (1, 0), (2, 0)]
j = [(0, 0), (0, 1), (0, 2), (1, 2)]

k = [(0, 0), (1, 0), (1, 1), (2, 1)]
l = [(0, 0), (0, 1), (1, -1), (1, 0)]

m = [(0, 0), (1, 0), (1, -1), (2, -1)]
n = [(0, 0), (0, 1), (1, 1), (1, 2)]

o = [(0, 0), (1, 0), (1, 1), (2, 0)]
p = [(0, 0), (1, -1), (1, 0), (1, 1)]
q = [(0, 0), (1, -1), (1, 0), (2, 0)]
r = [(0, 0), (0, 1), (0, 2), (1, 1)]

s = [(0, 0), (0, 1), (1, 0), (1, 1)]

tetrolookup = {"A": [a, "a"], "B": [b, "b"], "C": [c, "c"], "D": [d, "d"],
               "E": [e, "e"], "F": [f, "f"], "G": [g, "g"], "H": [h, "h"],
               "I": [i, "i"], "J": [j, "j"], "K": [k, "k"], "L": [l, "l"],
               "M": [m, "m"], "N": [n, "n"], "O": [o, "o"], "P": [p, "p"],
               "Q": [q, "q"], "R": [r, "r"], "S": [s, "s"]}

tetrolist = [tetrolookup[tetro] for tetro in tetrominolist]

# -----------------------------------------------------------------------------
print "Generating random tetromino partitions of a %i by %i grid..." % (rows, columns)

partitions = []
polyominolist = []

while len(partitions) < iterations:
    region = [list(row) for row in blank]
    polylist = []
    iteration = 1
    filled = False

    while True:
        cursor = False
        for rowindex in range(rows):
            for columnindex in range(columns):
                if region[rowindex][columnindex] == 0:
                    cursor = [rowindex, columnindex]
                    break
            if cursor:
                break
        
        if not cursor:
            filled = True
            break
            
        polyominos = list(tetrolist)
        random.shuffle(polyominos)
        polyplaced = False
        while len(polyominos) > 0:
            poly, polytext = polyominos.pop()
    
            if (cursor[0] + poly[0][0]) >= 0 and (cursor[1] + poly[0][1]) >= 0 and \
               (cursor[0] + poly[1][0]) >= 0 and (cursor[1] + poly[1][1]) >= 0 and \
               (cursor[0] + poly[2][0]) >= 0 and (cursor[1] + poly[2][1]) >= 0 and \
               (cursor[0] + poly[3][0]) >= 0 and (cursor[1] + poly[3][1]) >= 0 and \
               \
               (cursor[0] + poly[0][0]) < rows and (cursor[1] + poly[0][1]) < columns and \
               (cursor[0] + poly[1][0]) < rows and (cursor[1] + poly[1][1]) < columns and \
               (cursor[0] + poly[2][0]) < rows and (cursor[1] + poly[2][1]) < columns and \
               (cursor[0] + poly[3][0]) < rows and (cursor[1] + poly[3][1]) < columns and \
               \
               region[cursor[0] + poly[0][0]][cursor[1] + poly[0][1]] == 0 and \
               region[cursor[0] + poly[1][0]][cursor[1] + poly[1][1]] == 0 and \
               region[cursor[0] + poly[2][0]][cursor[1] + poly[2][1]] == 0 and \
               region[cursor[0] + poly[3][0]][cursor[1] + poly[3][1]] == 0:
    
                region[cursor[0] + poly[0][0]][cursor[1] + poly[0][1]] = iteration
                region[cursor[0] + poly[1][0]][cursor[1] + poly[1][1]] = iteration
                region[cursor[0] + poly[2][0]][cursor[1] + poly[2][1]] = iteration
                region[cursor[0] + poly[3][0]][cursor[1] + poly[3][1]] = iteration
                polyplaced = True
                break
        
        if polyplaced:
            iteration += 1
            polylist.append(polytext)
        else:
            break

    if filled:
        if region not in partitions:
            partitions.append(region)
            polyominolist.extend(polylist)
            if len(partitions) % 100 == 0:
                print len(partitions)

# -----------------------------------------------------------------------------
outputfile = open(pickledoutput, "wb")
pickle.dump(partitions, outputfile, -1)
outputfile.close()

outputfile = open(countoutput, "w")
for tetro in tetrominolist:
    outputfile.write("%s: %i\n" % (tetro.lower(), polyominolist.count(tetro.lower())))
outputfile.close()

# -----------------------------------------------------------------------------
hidlist = []
hidfile = open(idfile, "r")

for record in hidfile.readlines()[1:]:
    record = record.strip("\n").split("\t")
    hidlist.append(int(record[idfields.index("ID")]))

hidfile.close()
hidlist.sort()

centroids = {}
centroidfile = open(xyfile, "r")

for record in centroidfile.readlines()[1:]:
    record = record.strip("\n").split("\t")

    hid = int(record[xyfields.index("ID")])
    if hid in hidlist:
        centroids[hid] = Point(float(record[xyfields.index("X")]),
                               float(record[xyfields.index("Y")]))

centroidfile.close()

# -----------------------------------------------------------------------------
nx = len(hidlist) / columns
breaksx = range(0, len(hidlist), nx)
if len(hidlist) % columns == 0:
    breaksx.append(len(hidlist))
else:
    breaksx[-1] = len(hidlist)

origin = Point(0.0, 0.0)

# -----------------------------------------------------------------------------
print "Done. Applying tetromino partitions to input centroids..."

folds = {}
for hid in hidlist:
    folds[hid] = []

for iteration in range(iterations):
    rangle = random.uniform(anglelist[0], anglelist[1])
    xylist = []

    for hid in hidlist:
        rcentroid = affinity.rotate(centroids[hid], rangle, origin)
        xylist.append([rcentroid.x, rcentroid.y, hid])
        
    xylist.sort()
    for xindex in range(columns):
        ylist = []
        for hexagon in xylist[breaksx[xindex]:breaksx[xindex + 1]]:
            ylist.append([hexagon[1], hexagon[2]])

        ny = len(ylist) / rows
        breaksy = range(0, len(ylist), ny)
        if len(ylist) % rows == 0:
            breaksy.append(len(ylist))
        else:
            breaksy[-1] = len(ylist)

        ylist.sort()
        for yindex in range(rows):
            for hexagon in ylist[breaksy[yindex]:breaksy[yindex + 1]]:
                folds[hexagon[1]].append(partitions[iteration][yindex][xindex])

    if (iteration + 1) % 100 == 0:
        print iteration + 1

# -----------------------------------------------------------------------------
outputfile = open(foldsoutput, "w")

outputfile.write("ID")
for iteration in range(iterations):
    outputfile.write("\tI%04i" % (iteration + 1))

for hid in hidlist:
    outputfile.write("\n%i" % hid)
    for iteration in range(iterations):
        outputfile.write("\t%i" % folds[hid][iteration])

outputfile.close()
print "Done."
