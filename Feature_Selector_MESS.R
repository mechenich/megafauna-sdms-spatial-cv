arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))

species <- read.table(paste(s, "/", s, "_PA_Natural_O", n, ".txt", sep = ""),
                      header = TRUE, sep = "\t")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

terra <- read.table("../Ecosphere Analysis/Ecosphere Datasets/ISEA3H09_NE10M_V040100_Terra_Fractions.txt",
                    header = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
species <- merge(species, bio, by = "HID")
species <- merge(species, climdex, by = "HID")

terra <- merge(terra, bio, by = "HID")
terra <- merge(terra, climdex, by = "HID")
terra <- terra[terra$Terra_Fraction >= 0.5, ]
terra <- terra[terra$BIO01_Mean != -1000, ]

# -----------------------------------------------------------------------------
source("MESS_Functions.R")

features <- c(names(bio[, -1]), names(climdex[, -1]))

p.matrix <- ecdf.matrix(features, species, terra)
m.matrix <- mess.matrix(features, species, terra, p.matrix)
c.frame <- pearson.frame(features, species)

mess.global.frame <- messcor(features, m.matrix, c.frame)
mess.global.remove <- r.vector(mess.global.frame, 0.9)

f.set <- setdiff(features, mess.global.remove)
print(f.set)
save(f.set, file = paste(s, "/", s, "_Features_Natural_O", n, "_MESS9.rdata", sep = ""))

