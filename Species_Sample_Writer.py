import sys

s = sys.argv[1]
n = int(sys.argv[2])
cv = sys.argv[3]

columns = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 991, 992, 993, 994, 995, 996, 997, 998, 999, 1000]


inputpath = "%s/%s_Folds_%s10_Natural_O%02i.txt" % (s, s, cv, n)
outputpath = inputpath[:-4] + "_Sample.txt"

inputfile = open(inputpath, "r")
outputfile = open(outputpath, "w")

for record in inputfile.readlines():
    record = record.strip("\n").split("\t")
    
    outputtext = record[columns[0]]
    for column in columns[1:]:
        outputtext += "\t%s" % record[column]

    outputfile.write("%s\n" % outputtext)

inputfile.close()
outputfile.close()
