import pickle
import sys

o = sys.argv[1]
s = sys.argv[2]
e = float(sys.argv[3])
i = int(sys.argv[4])

# -----------------------------------------------------------------------------
dgg = {}

terrapath = "../Ecosphere Analysis/Ecosphere Datasets/" + \
            "ISEA3H09_NE10M_V040100_Terra_Fractions.txt"
terrafile = open(terrapath, "r")

for record in terrafile.readlines()[1:]:
    record = record.strip("\n").split("\t")
    dgg[int(record[0])] = [float(record[1])]

terrafile.close()

topopath = "../../Ecosphere/ISEA3H09/SRTM30PLUS_V11/" + \
           "ISEA3H09_SRTM30PLUS_V11_Elevation_Mean.txt"
topofile = open(topopath, "r")

for record in topofile.readlines()[1:]:
    record = record.strip("\n").split("\t")
    dgg[int(record[0])].append(float(record[1]))

topofile.close()

# -----------------------------------------------------------------------------
plist = []

rangepath = "../../Ecosphere/ISEA3H09/PHYLACINE_V010201/" + \
            "ISEA3H09_PHYLACINE_V010201_PresentNatural_%s_Centroid.txt"
rangefile = open(rangepath % o, "r")

header = rangefile.readline().upper().strip("\n").split("\t")
index = header.index("%s_CENTROID" % s.upper())

record = rangefile.readline()
while record:
    record = record.strip("\n").split("\t")
    hid = int(record[0])

    if record[index] == "1" and dgg[hid][0] >= 0.5 and dgg[hid][1] <= e:
        plist.append(hid)

    record = rangefile.readline()

rangefile.close()

# -----------------------------------------------------------------------------
alist = []

picklepath = "../Ecosphere Analysis/Spatial Statistics/Neighbors/" + \
             "Terra/Neighbors - ISEA3H09 - %02i.pkl"

for order in range(1, i + 1):
    picklefile = open(picklepath % order, "rb")
    neighbors = pickle.load(picklefile)
    picklefile.close()
    
    for phid in plist:
        for nhid in neighbors[phid]:
            
            if nhid not in plist and nhid not in alist and \
            dgg[nhid][1] <= e:
                alist.append(nhid)
    
    print order

# -----------------------------------------------------------------------------
outpath = "%s/%s_PA_Natural_O%02i.txt"

outfile = open(outpath % (s, s, i), "w")
outfile.write("HID\tPA")

for hid in range(1, 196833):
    if hid in plist:
        outfile.write("\n%i\t1" % hid)
    elif hid in alist:
        outfile.write("\n%i\t0" % hid)

outfile.close()

total = len(plist) + len(alist)
print "Presences: %i (%0.4f), Absences: %i, Total: %i" % (len(plist),
      (len(plist) / float(total)), len(alist), total)
