library(pROC)

setwd("/home/ad/home/m/mechenic/Projects/SDM Rewilding")

# -----------------------------------------------------------------------------
species <- read.table("Bison_Bonasus/Bison_Bonasus_PA_Natural_O20.txt",
                      header = TRUE, sep = "\t")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

species <- merge(species, bio, by = "HID")
species <- merge(species, climdex, by = "HID")

# -----------------------------------------------------------------------------
# folds <- read.table("Bison_Bonasus/Bison_Bonasus_Folds_R10_Natural_O20.txt",
folds <- read.table("Bison_Bonasus/Bison_Bonasus_Folds_S10_Natural_O20.txt",
                    header = TRUE, sep = "\t")[, -1]
# 0.9
vifcor <- c("BIO02_Mean", "BIO03_Mean", "BIO05_Mean", "BIO07_Mean", "BIO08_Mean",
            "BIO09_Mean", "BIO12_Mean", "BIO13_Mean", "BIO14_Mean", "BIO15_Mean",
            "BIO18_Mean", "BIO19_Mean", "CDD_IDW1N10", "CSDI_IDW1N10", "CWD_IDW1N10",
            "DTR_IDW1N10", "GSL_IDW1N10", "R10MM_IDW1N10", "R1MM_IDW1N10", "R99P_IDW1N10",
            "SDII_IDW1N10", "TN10P_IDW1N10", "TN90P_IDW1N10", "TNX_IDW1N10", "TR_IDW1N10",
            "TX10P_IDW1N10", "TX90P_IDW1N10", "TXX_IDW1N10", "WSDI_IDW1N10")

messcor <- c() 

# -----------------------------------------------------------------------------
for (i in 184:200) {
  features.in <- vifcor   # messcor
  features.out <- c()

  t.out <- c()
  v.out <- c()

  while (length(features.in) > 0) {
    t.features <- c()
    v.features <- c()

    for (v in features.in) {
      t.obs <- c()
      v.obs <- c()
      t.pre <- c()
      v.pre <- c()

      for (f in 1:10) {
        d.train <- species[folds[, i] != f, c("PA", v, features.out)]
        d.valid <- species[folds[, i] == f, c("PA", v, features.out)]

        m.fit <- glm(PA ~ ., data = d.train, family = "binomial")
        m.pre <- predict(m.fit, newdata = d.valid, type = "response")

        t.obs = c(t.obs, d.train$PA)
        v.obs = c(v.obs, d.valid$PA)
        t.pre = c(t.pre, m.fit$fitted.values)
        v.pre = c(v.pre, m.pre)
      }

      t.features <- c(t.features, auc(t.obs, t.pre, quiet = TRUE))
      v.features <- c(v.features, auc(v.obs, v.pre, quiet = TRUE))
    }

    index <- which.max(v.features)
    features.out <- c(features.out, features.in[index])
    features.in <- features.in[-1 * index]

    t.out <- c(t.out, t.features[index])
    v.out <- c(v.out, v.features[index])

    print(paste("Random", sprintf("%i", i),
                "- Feature", sprintf("%i", length(features.out)),
                "-", features.out[length(features.out)],
                "- AUC", sprintf("%0.4f", v.features[index])))
  }

  results <- data.frame(Feature.Index = 1:length(features.out),
                        Iteration = i,
                        Feature = features.out,
                        AUC.Train = t.out,
                        AUC.Valid = v.out,
                        stringsAsFactors = FALSE)
  
# save(results, file = paste("Bison_Bonasus_Natural_O20_R10_VIF9_GLM_AUC/",
#                            "Bison_Bonasus_Natural_O20_R10_VIF9_GLM_AUC_I",
save(results, file = paste("Bison_Bonasus_Natural_O20_S10_VIF9_GLM_AUC/",
                           "Bison_Bonasus_Natural_O20_S10_VIF9_GLM_AUC_I",
                           sprintf("%03i", i), ".rdata", sep = ""))
}
