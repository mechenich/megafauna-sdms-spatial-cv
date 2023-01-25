library(usdm)

arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))

species <- read.table(paste(s, "/", s, "_PA_Natural_O", n, ".txt", sep = ""),
                      header = TRUE, sep = "\t")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

species <- merge(species, bio, by = "HID")
species <- merge(species, climdex, by = "HID")

# -----------------------------------------------------------------------------
f.set <- vifcor(species[, 3:48], 0.9)@results$Variables
save(f.set, file = paste(s, "/", s, "_Features_Natural_O", n, "_VIF9.rdata", sep = ""))

# -----------------------------------------------------------------------------
r.matrix <- cor(species[, 3:48], method = "pearson")
r.matrix.text <- matrix(sprintf("%0.6f", r.matrix),
                        nrow = nrow(r.matrix), ncol = ncol(r.matrix))

row.names(r.matrix.text) <- row.names(r.matrix)
colnames(r.matrix.text) <- colnames(r.matrix)

write.table(r.matrix.text,
            file = paste(s, "/", s, "_R_Natural_O", n, ".txt", sep = ""),
            quote = FALSE, sep = "\t")
