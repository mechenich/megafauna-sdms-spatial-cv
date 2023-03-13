library(pROC)
setwd("/home/ad/home/m/mechenic/Projects/megafauna-sdms-spatial-cv")

neighbors <- c("30", "20", "30", "06", "15", "20", "20", "11", "06", "04")

i.start <- data.frame(id = c("VIF9-R10", "VIF9-S10", "MESS9-R10", "MESS9-S10"),
                      index = c(600, 700, 800, 900))

specieslist <- c("Alces_Alces", "Bison_Bonasus", "Camelus_Dromedarius",
                 "Ceratotherium_Simum", "Cervus_Elaphus", "Diceros_Bicornis",
                 "Elephas_Maximus", "Equus_Africanus", "Ovibos_Moschatus",
                 "Rhinoceros_Unicornis")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

# -----------------------------------------------------------------------------
for (s in specieslist) {
    species <- read.table(paste(s, "/", s, "_PA_Natural_O",
                                neighbors[match(s, specieslist)], ".txt", sep = ""),
                          header = TRUE, sep = "\t")

    species <- merge(species, bio, by = "HID")
    species <- merge(species, climdex, by = "HID")

    folds <- read.table(paste(s, "/", s, "_Folds_S10_Natural_O",
                              neighbors[match(s, specieslist)], ".txt", sep = ""),
                    header = TRUE, sep = "\t")[, -1]

    # -------------------------------------------------------------------------
    for (p in c("VIF9", "MESS9")) {
        for (cv in c("R10", "S10")) {
            load(paste(s, "/", s, "_Features_Natural_O",
                       neighbors[match(s, specieslist)], "_", cv, "_", p,
                       "_GLM.rdata", sep = ""))
          
            # -----------------------------------------------------------------
            e.id <- paste(p, "-", cv, sep = "")
            e.start <- i.start$index[i.start$id == e.id]
            e.auc <- c()

            for (i in e.start:(e.start + 100)) {
                v.obs <- c()
                v.pre <- c()
          
                for (f in 1:10) {
                    d.train <- species[folds[, i] != f, c("PA", f.set$Feature)]
                    d.valid <- species[folds[, i] == f, c("PA", f.set$Feature)]

                    m.fit <- glm(PA ~ ., data = d.train, family = "binomial")
                    m.pre <- predict(m.fit, newdata = d.valid, type = "response")

                    v.obs = c(v.obs, d.valid$PA)
                    v.pre = c(v.pre, m.pre)
                }
                e.auc <- c(e.auc, auc(v.obs, v.pre, quiet = TRUE))
            }
            save(e.auc, file = paste(s, "/", s, "_AUC_Natural_O",
                                     neighbors[match(s, specieslist)], "_", cv,
                                     "_", p, "_GLM.rdata", sep = ""))

            print(paste(s, "-", e.id, "-", sprintf("%.2f", median(e.auc)), sep = ""))
        }
    }
}