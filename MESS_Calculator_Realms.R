arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))
p <- toupper(arguments[3])
m <- toupper(arguments[4])

# -----------------------------------------------------------------------------
species <- read.table(paste(s, "/", s, "_PA_Natural_O", n, ".txt", sep = ""),
                      header = TRUE, sep = "\t")

terra <- read.table("../Ecosphere Analysis/Ecosphere Datasets/ISEA3H09_NE10M_V040100_Terra_Fractions.txt",
                    header = TRUE, sep = "\t")

realm <- read.table("../../Ecosphere/ISEA3H09/WWFTE_V02/ISEA3H09_WWFTE_V02_Realm_Mode.txt",
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE)

realm$Realm_Mode[is.na(realm$Realm_Mode)]  <- "Nearctic"
realm$Realm_Mode[realm$Realm_Mode == "PA"] <- "Palearctic"
realm$Realm_Mode[realm$Realm_Mode == "OC"] <- "Oceania"
realm$Realm_Mode[realm$Realm_Mode == "NT"] <- "Neotropic"
realm$Realm_Mode[realm$Realm_Mode == "AT"] <- "Afrotropic"
realm$Realm_Mode[realm$Realm_Mode == "IM"] <- "Indo-Malay"
realm$Realm_Mode[realm$Realm_Mode == "AA"] <- "Australasia"
realm$Realm_Mode[realm$Realm_Mode == "AN"] <- "Antarctic"
realm$Realm_Mode[realm$Realm_Mode == "-1"] <- "None"

realms <- unique(realm$Realm_Mode)
realms <- realms[order(realms)]

# -----------------------------------------------------------------------------
bio.1950 <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                       header = TRUE, sep = "\t")
bio.2061 <- read.table("WC30AS_V14_CMIP5_CCSM4_RCP85_2070_BIO/ISEA3H09_WC30AS_V14_CMIP5_CCSM4_RCP85_2070_BIO.txt",
                       header = TRUE, sep = "\t")

climdex.1950 <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                           header = TRUE, sep = "\t")
climdex.2061 <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y2061_Y2080_ETCCDI_IDW1N10.txt",
                           header = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
trimmer <- function(text) {strsplit(text, "_")[[1]][1]}

features <- c(sapply(names(bio.1950[, -1]), trimmer),
              sapply(names(climdex.1950[, -1]), trimmer))

species <- merge(species, bio.1950, by = "HID")
species <- merge(species, climdex.1950, by = "HID")
species <- species[, -2]
names(species) <- c("HID", features)

terra <- terra[terra$Terra_Fraction >= 0.5, ]
terra <- merge(terra, realm, by = "HID")

terra.1950 <- merge(terra, bio.1950, by = "HID")
terra.1950 <- merge(terra.1950, climdex.1950, by = "HID")
terra.1950 <- terra.1950[terra.1950$BIO01_Mean != -1000, ]
terra.1950 <- terra.1950[, -2]
names(terra.1950) <- c("HID", "Realm", features)

terra.2061 <- merge(terra, bio.2061, by = "HID")
terra.2061 <- merge(terra.2061, climdex.2061, by = "HID")
terra.2061 <- terra.2061[terra.2061$BIO01_Mean != -1000, ]
terra.2061 <- terra.2061[, -2]
names(terra.2061) <- c("HID", "Realm", features)

# -----------------------------------------------------------------------------
f.sets <- data.frame(ID = integer(), Feature = character(), Feature.Short = character(),
                     Frequency = double(), CV = character(), stringsAsFactors = FALSE)
for (cv in c("R", "S")) {
  load(paste(s, "/", s, "_Features_Natural_O", n, "_", cv, "10_", p, "_", m, ".rdata", sep = ""))
  f.set$CV <- cv
  f.sets <- rbind(f.sets, f.set)
}

feature.set <- unique(f.sets$Feature.Short)
r.set <- feature.set %in% f.sets[f.sets$CV == "R", "Feature.Short"]
s.set <- feature.set %in% f.sets[f.sets$CV == "S", "Feature.Short"]

# -----------------------------------------------------------------------------
source("MESS_Functions.R")

realm.si.1950 <- matrix(nrow = length(realms), ncol = length(feature.set))
realm.sip.1950 <- matrix(nrow = length(realms), ncol = length(feature.set))
realm.si.2061 <- matrix(nrow = length(realms), ncol = length(feature.set))
realm.sip.2061 <- matrix(nrow = length(realms), ncol = length(feature.set))

realm.m.r.1950 <- matrix(nrow = 2, ncol = length(realms))
realm.m.s.1950 <- matrix(nrow = 2, ncol = length(realms))
realm.m.r.2061 <- matrix(nrow = 2, ncol = length(realms))
realm.m.s.2061 <- matrix(nrow = 2, ncol = length(realms))

for (r in realms) {
    p.matrix.1950 <- ecdf.matrix(feature.set,
                                 species[, feature.set],
                                 terra.1950[terra.1950$Realm == r, feature.set])
    p.matrix.2061 <- ecdf.matrix(feature.set,
                                 species[, feature.set],
                                 terra.2061[terra.2061$Realm == r, feature.set])

    m.matrix.1950 <- mess.matrix(feature.set,
                                 species[, feature.set],
                                 terra.1950[terra.1950$Realm == r, feature.set],
                                 p.matrix.1950)
    m.matrix.2061 <- mess.matrix(feature.set,
                                 species[, feature.set],
                                 terra.2061[terra.2061$Realm == r, feature.set],
                                 p.matrix.2061)
    
    realm.si.1950[match(r, realms), ] <- apply(m.matrix.1950, 2, mean)
    realm.sip.1950[match(r, realms), ] <- apply(m.matrix.1950 > 0, 2, sum) / nrow(m.matrix.1950)
    
    realm.m.r.1950[1, match(r, realms)] <- mean(apply(m.matrix.1950[, r.set], 1, min))
    realm.m.r.1950[2, match(r, realms)] <- sum(apply(m.matrix.1950[, r.set], 1, min) > 0) / nrow(m.matrix.1950)
    realm.m.s.1950[1, match(r, realms)] <- mean(apply(m.matrix.1950[, s.set], 1, min))
    realm.m.s.1950[2, match(r, realms)] <- sum(apply(m.matrix.1950[, s.set], 1, min) > 0) / nrow(m.matrix.1950)

    realm.si.2061[match(r, realms), ] <- apply(m.matrix.2061, 2, mean)
    realm.sip.2061[match(r, realms), ] <- apply(m.matrix.2061 > 0, 2, sum) / nrow(m.matrix.2061)

    realm.m.r.2061[1, match(r, realms)] <- mean(apply(m.matrix.2061[, r.set], 1, min))
    realm.m.r.2061[2, match(r, realms)] <- sum(apply(m.matrix.2061[, r.set], 1, min) > 0) / nrow(m.matrix.2061)
    realm.m.s.2061[1, match(r, realms)] <- mean(apply(m.matrix.2061[, s.set], 1, min))
    realm.m.s.2061[2, match(r, realms)] <- sum(apply(m.matrix.2061[, s.set], 1, min) > 0) / nrow(m.matrix.2061)
}

# -----------------------------------------------------------------------------
p.matrix.1950 <- ecdf.matrix(feature.set, species[, feature.set], terra.1950[, feature.set])
p.matrix.2061 <- ecdf.matrix(feature.set, species[, feature.set], terra.2061[, feature.set])
  
m.matrix.1950 <- mess.matrix(feature.set, species[, feature.set],
                             terra.1950[, feature.set], p.matrix.1950)
m.matrix.2061 <- mess.matrix(feature.set, species[, feature.set],
                             terra.2061[, feature.set], p.matrix.2061)
  
global.si.1950 <- apply(m.matrix.1950, 2, mean)
global.sip.1950 <- apply(m.matrix.1950 > 0, 2, sum) / nrow(m.matrix.1950)
global.m.r.1950 <- mean(apply(m.matrix.1950[, r.set], 1, min))
global.mp.r.1950 <- sum(apply(m.matrix.1950[, r.set], 1, min) > 0) / nrow(m.matrix.1950)
global.m.s.1950 <- mean(apply(m.matrix.1950[, s.set], 1, min))
global.mp.s.1950 <- sum(apply(m.matrix.1950[, s.set], 1, min) > 0) / nrow(m.matrix.1950)
  
global.si.2061 <- apply(m.matrix.2061, 2, mean)
global.sip.2061 <- apply(m.matrix.2061 > 0, 2, sum) / nrow(m.matrix.2061)
global.m.r.2061 <- mean(apply(m.matrix.2061[, r.set], 1, min))
global.mp.r.2061 <- sum(apply(m.matrix.2061[, r.set], 1, min) > 0) / nrow(m.matrix.2061)
global.m.s.2061 <- mean(apply(m.matrix.2061[, s.set], 1, min))
global.mp.s.2061 <- sum(apply(m.matrix.2061[, s.set], 1, min) > 0) / nrow(m.matrix.2061)

# -----------------------------------------------------------------------------
si.r.1950 <- c(mean(global.si.1950[r.set]), apply(realm.si.1950[, r.set], 1, mean))
sip.r.1950 <- c(mean(global.sip.1950[r.set]), apply(realm.sip.1950[, r.set], 1, mean))
si.s.1950 <- c(mean(global.si.1950[s.set]), apply(realm.si.1950[, s.set], 1, mean))
sip.s.1950 <- c(mean(global.sip.1950[s.set]), apply(realm.sip.1950[, s.set], 1, mean))

si.1950 <- rbind(si.r.1950, sip.r.1950, si.s.1950, sip.s.1950)
colnames(si.1950) <- c("Global", realms)
row.names(si.1950) <- c("R-SI", "R-%", "S-SI", "S-%")

write.table(si.1950,
            paste(s, "/", s, "_Realms_Y1950_SI_Natural_O", n, "_", p, "_", m, ".txt",
                  sep = ""), sep = "\t", quote = FALSE)

si.r.2061 <- c(mean(global.si.2061[r.set]), apply(realm.si.2061[, r.set], 1, mean))
sip.r.2061 <- c(mean(global.sip.2061[r.set]), apply(realm.sip.2061[, r.set], 1, mean))
si.s.2061 <- c(mean(global.si.2061[s.set]), apply(realm.si.2061[, s.set], 1, mean))
sip.s.2061 <- c(mean(global.sip.2061[s.set]), apply(realm.sip.2061[, s.set], 1, mean))

si.2061 <- rbind(si.r.2061, sip.r.2061, si.s.2061, sip.s.2061)
colnames(si.2061) <- c("Global", realms)
row.names(si.2061) <- c("R-SI", "R-%", "S-SI", "S-%")

write.table(si.2061,
            paste(s, "/", s, "_Realms_Y2061_SI_Natural_O", n, "_", p, "_", m, ".txt",
                  sep = ""), sep = "\t", quote = FALSE)

m.r.1950 <- cbind(c(global.m.r.1950, global.mp.r.1950), realm.m.r.1950)
m.s.1950 <- cbind(c(global.m.s.1950, global.mp.s.1950), realm.m.s.1950)
m.1950 <- rbind(m.r.1950, m.s.1950)

colnames(m.1950) <- c("Global", realms)
row.names(m.1950) <- c("R-MESS", "R-%", "S-MESS", "S-%")

write.table(m.1950,
            paste(s, "/", s, "_Realms_Y1950_MESS_Natural_O", n, "_", p, "_", m, ".txt",
            sep = ""), sep = "\t", quote = FALSE)

m.r.2061 <- cbind(c(global.m.r.2061, global.mp.r.2061), realm.m.r.2061)
m.s.2061 <- cbind(c(global.m.s.2061, global.mp.s.2061), realm.m.s.2061)
m.2061 <- rbind(m.r.2061, m.s.2061)

colnames(m.2061) <- c("Global", realms)
row.names(m.2061) <- c("R-MESS", "R-%", "S-MESS", "S-%")

write.table(m.2061,
            paste(s, "/", s, "_Realms_Y2061_MESS_Natural_O", n, "_", p, "_", m, ".txt",
            sep = ""), sep = "\t", quote = FALSE)
