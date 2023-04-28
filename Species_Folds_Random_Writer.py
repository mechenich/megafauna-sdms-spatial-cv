import random
import sys

s = sys.argv[1]
n = int(sys.argv[2])
i = int(sys.argv[3])

r = 0
for index in range(len(s)):
    r += ord(s[index])
print r

# -----------------------------------------------------------------------------
hidlist = []

speciesfile = open("%s/%s_PA_Natural_O%02i.txt" % (s, s, n), "r")

for record in speciesfile.readlines()[1:]:
    record = record.strip("\n").split("\t")
    hidlist.append(int(record[0]))

speciesfile.close()

# -----------------------------------------------------------------------------
k = 10
interval = len(hidlist) / k

bins = range(0, len(hidlist), interval)
if len(hidlist) % k == 0:
    bins.append(len(hidlist))
else:
    bins[-1] = len(hidlist)

folds = {}
for hid in hidlist:
    folds[hid] = []

# -----------------------------------------------------------------------------
random.seed(r)

for iteration in range(i):
    random.shuffle(hidlist)
    
    for index in range(k):
        for hid in hidlist[bins[index]:bins[index + 1]]:
            folds[hid].append(index + 1)

# -----------------------------------------------------------------------------
hidlist.sort()

outfile = open("%s/%s_Folds_R10_Natural_O%02i.txt" % (s, s, n), "w")

outfile.write("HID")
for iteration in range(i):
    outfile.write("\tI%03i" % (iteration + 1))

for hid in hidlist:
    outfile.write("\n%i" % hid)
    for iteration in range(i):
        outfile.write("\t%i" % folds[hid][iteration])

outfile.close()
