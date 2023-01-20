library(usdm)

setwd("/home/ad/home/m/mechenic/Projects/SDM Rewilding")

bison <- read.table("Bison_Bonasus/Bison_Bonasus_PA_Natural_O20.txt",
                    header = TRUE, sep = "\t")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

bison <- merge(bison, bio, by = "HID")
bison <- merge(bison, climdex, by = "HID")

# -----------------------------------------------------------------------------
x <- vifcor(bison[, 3:48], 0.9)
x@results$Variables

# -----------------------------------------------------------------------------
r.matrix <- cor(bison[, 3:48], method = "pearson")
r.matrix.text <- matrix(sprintf("%0.6f", r.matrix),
                        nrow = nrow(r.matrix), ncol = ncol(r.matrix))

row.names(r.matrix.text) <- row.names(r.matrix)
colnames(r.matrix.text) <- colnames(r.matrix)

write.table(r.matrix.text,
            file = "Bison_Bonasus/Bison_Bonasus_R_Natural_O20.txt",
            quote = FALSE, sep = "\t")
