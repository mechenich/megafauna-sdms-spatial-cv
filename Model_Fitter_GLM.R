library(pROC)

arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))
cv <- toupper(arguments[3])
istart <- as.numeric(arguments[4])
istop <- as.numeric(arguments[5])
p <- toupper(arguments[6])

# -----------------------------------------------------------------------------
species <- read.table(paste(s, "/", s, "_PA_Natural_O", n, ".txt", sep = ""),
                      header = TRUE, sep = "\t")

bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

species <- merge(species, bio, by = "HID")
species <- merge(species, climdex, by = "HID")

# -----------------------------------------------------------------------------
folds <- read.table(paste(s, "/", s, "_Folds_", cv, "10_Natural_O", n, ".txt", sep = ""),
                    header = TRUE, sep = "\t")[, -1]

load(paste(s, "/", s, "_Features_Natural_O", n, "_", p, ".rdata", sep = ""))

# -----------------------------------------------------------------------------
for (i in istart:istop) {
  features.in <- f.set
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

    print(paste("Iteration", sprintf("%i", i),
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
  
save(results, file = paste(s, "_Natural_O", n, "_", cv, "10_", p, "_GLM_AUC/",
                           s, "_Natural_O", n, "_", cv, "10_", p, "_GLM_AUC_I",
                           sprintf("%03i", i), ".rdata", sep = ""))
}
