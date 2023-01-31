arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))
p <- toupper(arguments[3])
m <- toupper(arguments[4])

# -----------------------------------------------------------------------------
bio <- read.table("WC30AS_V14_BIO/ISEA3H09_WC30AS_V14_BIO.txt",
                  header = TRUE, sep = "\t")

climdex <- read.table("CCSM4_ETCCDI/ISEA3H09_CCSM4_Y1950_Y2000_ETCCDI_IDW1N10.txt",
                      header = TRUE, sep = "\t")

folds <- read.table(paste(s, "/", s, "_Folds_S10_Natural_O", n, ".txt", sep = ""),
                    header = TRUE, sep = "\t")

trimmer <- function(text) {strsplit(text, "_")[[1]][1]}

features <- c(sapply(names(bio[, -1]), trimmer),
              sapply(names(climdex[, -1]), trimmer))

training <- data.frame(HID = folds[, 1])
training <- merge(training, bio, by = "HID")
training <- merge(training, climdex, by = "HID")
names(training) <- c("HID", features)

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
extract.i <- function(text) {
    as.numeric(strsplit(strsplit(text, "_I")[[1]][2], "\\.")[[1]][1])
}

iterations <- c(sapply(list.files(paste(s, "_Natural_O", n, "_S10_", p, "_", m,
                                        "_AUC/", sep = "")), extract.i))

# -----------------------------------------------------------------------------
source("MESS_Functions.R")

si.mean <- matrix(nrow = length(iterations), ncol = length(feature.set))
si.min <- matrix(nrow = length(iterations), ncol = length(feature.set))
si.p.mean <- matrix(nrow = length(iterations), ncol = length(feature.set))
si.p.min <- matrix(nrow = length(iterations), ncol = length(feature.set))
r.mess <- matrix(nrow = length(iterations), ncol = 4)
s.mess <- matrix(nrow = length(iterations), ncol = 4)

for (i in iterations) {
    folds.i <- folds[, i + 1]

    si.i <- matrix(nrow = 10, ncol = length(feature.set))
    si.p.i <- matrix(nrow = 10, ncol = length(feature.set))

    r.mess.i <- numeric(10)
    r.mess.p.i <- numeric(10)
    s.mess.i <- numeric(10)
    s.mess.p.i <- numeric(10)

    for (k in 1:10) {
        training.k <- training[folds.i %in% c(1:10)[-k], feature.set]
        valid.k <- training[folds.i == k, feature.set]
        p.matrix <- ecdf.matrix(feature.set, training.k, valid.k)
        m.matrix <- mess.matrix(feature.set, training.k, valid.k, p.matrix)
        
        for (f in 1:length(feature.set)) {
            si.i[k, f] <- mean(m.matrix[, f])
            si.p.i[k, f] <- sum(m.matrix[, f] > 0) / nrow(m.matrix)
        }

        r.mess.v <- apply(m.matrix[, r.set], 1, min)
        s.mess.v <- apply(m.matrix[, s.set], 1, min)
        r.mess.i[k] <- mean(r.mess.v)
        s.mess.i[k] <- mean(s.mess.v)
        r.mess.p.i[k] <- sum(r.mess.v > 0) / nrow(m.matrix)
        s.mess.p.i[k] <- sum(s.mess.v > 0) / nrow(m.matrix)
    }

    si.mean[match(i, iterations), ] <- apply(si.i, 2, mean) # Mean of folds.
    si.min[match(i, iterations), ] <- apply(si.i, 2, min)   # Minimum of folds.
    si.p.mean[match(i, iterations), ] <- apply(si.p.i, 2, mean)
    si.p.min[match(i, iterations), ] <- apply(si.p.i, 2, min)

    r.mess[match(i, iterations), ] <- c(mean(r.mess.i), mean(r.mess.p.i), min(r.mess.i), min(r.mess.p.i))
    s.mess[match(i, iterations), ] <- c(mean(s.mess.i), mean(s.mess.p.i), min(s.mess.i), min(s.mess.p.i))

    print(sprintf("%i", i))
}

# -----------------------------------------------------------------------------
si.results <- data.frame(Feature = feature.set,
                         R = r.set,
                         S = s.set,
                         R.Frequency = 0.0,
                         S.Frequency = 0.0,
                         SI.Mean = apply(si.mean, 2, mean),
                         SI.P.Mean = apply(si.p.mean, 2, mean),
                         SI.Min = apply(si.min, 2, mean),
                         SI.P.Min = apply(si.p.min, 2, mean))

for (f in si.results[si.results$R, "Feature"]) {
    si.results[si.results$Feature == f, "R.Frequency"] <- f.sets[f.sets$CV == "R" &
                                                                 f.sets$Feature.Short == f,
                                                                 "Frequency"]
}

for (f in si.results[si.results$S, "Feature"]) {
    si.results[si.results$Feature == f, "S.Frequency"] <- f.sets[f.sets$CV == "S" &
                                                                 f.sets$Feature.Short == f,
                                                                 "Frequency"]
}

r.means <- c("R.Mean", FALSE, FALSE, 0.0, 0.0, 
             mean(si.results[si.results$R, "SI.Mean"]),
             mean(si.results[si.results$R, "SI.P.Mean"]),
             mean(si.results[si.results$R, "SI.Min"]),
             mean(si.results[si.results$R, "SI.P.Min"]))

s.means <- c("S.Mean", FALSE, FALSE, 0.0, 0.0,
             mean(si.results[si.results$S, "SI.Mean"]),
             mean(si.results[si.results$S, "SI.P.Mean"]),
             mean(si.results[si.results$S, "SI.Min"]),
             mean(si.results[si.results$S, "SI.P.Min"]))

write.table(rbind(si.results, r.means, s.means),
            paste(s, "/", s, "_Folds_SI_Natural_O", n, "_", p, "_", m, ".txt",
                  sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)

m.results <- data.frame(CV = c("R", "S"),
                        M.Mean = c(mean(r.mess[, 1]), mean(s.mess[, 1])),
                        M.P.Mean = c(mean(r.mess[, 2]), mean(s.mess[, 2])),
                        M.Min = c(mean(r.mess[, 3]), mean(s.mess[, 3])),
                        M.P.Min = c(mean(r.mess[, 4]), mean(s.mess[, 4])))

write.table(m.results,
            paste(s, "/", s, "_Folds_MESS_Natural_O", n, "_", p, "_", m, ".txt",
                  sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
