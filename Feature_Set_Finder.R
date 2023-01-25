arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))
p <- toupper(arguments[3])
m <- toupper(arguments[4])
f <- as.numeric(arguments[5])

for (cv in c("R", "S")) {
    r.directory <- paste(s, "_Natural_O", n, "_", cv, "10_", p, "_", m, "_AUC/",
                         sep = "")
    r.files <- list.files(r.directory)

    r.frame <- data.frame(Feature.Index = integer(), Iteration = integer(),
                          Feature = character(), AUC.Train = double(),
                          AUC.Valid = double(), stringsAsFactors = FALSE)

    for (r in r.files) {
        load(paste(r.directory, r, sep = ""))
        r.frame <- rbind(r.frame, results)
    }

    r.frame$FF <- factor(r.frame$Feature)
    min.i <- min(r.frame$Iteration)
    max.i <- max(r.frame$Iteration)
    range.i <- length(min.i:max.i)
    offset.i <- 1 - min.i

    feature.matrix <- matrix(nrow = range.i, ncol = f)
    for (iteration in min.i:max.i) {
        for (feature in 1:f) {
            feature.matrix[iteration + offset.i, feature] <- r.frame$FF[r.frame$Iteration == iteration &
                                                                        r.frame$Feature.Index == feature]
        }
    }

    counter <- function(id) {
        sum(feature.matrix == id) / range.i
    }

    shortener <- function(text) {
        strsplit(text, "_")[[1]][1]
    }

    f.set <- data.frame(ID = 1:max(r.frame$Feature.Index))
    f.set$Feature <- levels(r.frame$FF)
    f.set$Feature.Short <- sapply(f.set$Feature, shortener)
    f.set$Frequency <- sapply(f.set$ID, counter)

    f.set <- f.set[order(f.set$Frequency, decreasing = TRUE), ]
    f.set <- f.set[1:f, ]

    print(f.set[, 3:4])
    save(f.set, file = paste(s, "/", s, "_Features_Natural_O", n, "_", cv, "10_",
                             p, "_", m, ".rdata", sep = ""))
}