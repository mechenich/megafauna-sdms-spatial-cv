z <- "/home/ad/home/m/mechenic/Projects/megafauna-sdms-spatial-cv"
setwd(z)

neighbors <- c("30", "20", "30", "06", "15", "20", "20", "11", "06", "04")

specieslist <- c("Alces_Alces", "Bison_Bonasus", "Camelus_Dromedarius",
                 "Ceratotherium_Simum", "Cervus_Elaphus", "Diceros_Bicornis",
                 "Elephas_Maximus", "Equus_Africanus", "Ovibos_Moschatus",
                 "Rhinoceros_Unicornis")

# ---------------------------------------------------------------------------------------
bio.1950 <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                       header = TRUE, sep = "\t")
bio.1950 <- bio.1950[bio.1950$BIO01_Mean > -1000, ]

bio.2061 <- read.table("WC30AS_V14_CMIP5_CCSM4_RCP85_2070_BIO/ISEA3H09_WC30AS_V14_CMIP5_CCSM4_RCP85_2070_BIO.txt",
                       header = TRUE, sep = "\t")
bio.2061 <- bio.2061[bio.2061$BIO01_Mean > -1000, ]

climdex.1950 <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                           header = TRUE, sep = "\t")

climdex.2061 <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y2061_Y2080_ETCCDI_IDW1N10.txt",
                           header = TRUE, sep = "\t")

climate.1950 <- merge(bio.1950, climdex.1950, by = "HID")
climate.2061 <- merge(bio.2061, climdex.2061, by = "HID")
sum(climate.1950$HID == climate.2061$HID)

# ---------------------------------------------------------------------------------------
for (s in specieslist) {
    s.predictions <- sprintf("%i", climate.1950$HID)

    species <- read.table(paste(s, "/", s, "_PA_Natural_O",
                                neighbors[match(s, specieslist)], ".txt", sep = ""),
                          header = TRUE, sep = "\t")

    species <- merge(species, climate.1950, by = "HID")
    sum(is.na(species$BIO01_Mean))
  
    for (p in c("VIF9", "MESS9")) {
        for (cv in c("R10", "S10")) {
            load(paste(s, "/", s, "_Features_Natural_O",
                       neighbors[match(s, specieslist)], "_", cv, "_", p,
                       "_GLM.rdata", sep = ""))

            m.fit <- glm(PA ~ ., data = species[, c("PA", f.set$Feature)], family = "binomial")

            m.1950 <- predict(m.fit, newdata = climate.1950, type = "response")
            m.2061 <- predict(m.fit, newdata = climate.2061, type = "response")
            s.predictions <- cbind(s.predictions, sprintf("%0.6f", m.1950))
            s.predictions <- cbind(s.predictions, sprintf("%0.6f", m.2061))
        }
    }

    colnames(s.predictions) <- c("HID", "VIF9-R10-1950", "VIF9-R10-2061",
                                        "VIF9-S10-1950", "VIF9-S10-2061",
                                        "MESS9-R10-1950", "MESS9-R10-2061",
                                        "MESS9-S10-1950", "MESS9-S10-2061")
    write.table(as.data.frame(s.predictions),
                file = paste(s, "/", s, "_Predictions_Natural_O", neighbors[match(s, specieslist)],
                             "_GLM.txt", sep = ""),
                quote = FALSE, sep = "\t", row.names = FALSE)
}
