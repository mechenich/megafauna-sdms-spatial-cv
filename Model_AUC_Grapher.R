library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
s <- arguments[1]
n <- sprintf("%02i", as.numeric(arguments[2]))
pool <- toupper(arguments[3])
m <- toupper(arguments[4])

results.frame <- data.frame(Feature.Index = integer(), Iteration = integer(),
                            Feature = character(), AUC.Train = double(),
                            AUC.Valid = double(), CV = character(), stringsAsFactors = FALSE)

for (cv in c("R", "S")) {
    r.directory <- paste(s, "_Natural_O", n, "_", cv, "10_", pool, "_", m, "_AUC/",
                         sep = "")
    r.files <- list.files(r.directory)

    for (r in r.files) {
        load(paste(r.directory, r, sep = ""))
        results$CV <- cv
        results.frame <- rbind(results.frame, results)
    }
}

feature.count <- max(results.frame$Feature.Index)

auc.random <- matrix(nrow = feature.count, ncol = 5)
auc.spatial <- matrix(nrow = feature.count, ncol = 5)

for (feature in 1:feature.count) {
    feature.auc <- results.frame[results.frame$Feature.Index == feature, c("AUC.Valid", "CV")]

    auc.random[feature, ] <- quantile(feature.auc[feature.auc$CV == "R", "AUC.Valid"],
                                      probs = seq(0, 1, 0.25))
    auc.spatial[feature, ] <- quantile(feature.auc[feature.auc$CV == "S", "AUC.Valid"],
                                       probs = seq(0, 1, 0.25))
}

feature.random <- which.max(auc.random[, 3])
feature.spatial <- which.max(auc.spatial[, 3])

p <- ggplot()
p <- p + geom_ribbon(aes(x = 1:feature.count, ymin = auc.random[, 1], ymax = auc.random[, 5]),
                     fill = "#fcbba180")
p <- p + geom_ribbon(aes(x = 1:feature.count, ymin = auc.random[, 2], ymax = auc.random[, 4]),
                     fill = "#fb6a4a80")

p <- p + geom_ribbon(aes(x = 1:feature.count, ymin = auc.spatial[, 1], ymax = auc.spatial[, 5]),
                     fill = "#c6dbef80")
p <- p + geom_ribbon(aes(x = 1:feature.count, ymin = auc.spatial[, 2], ymax = auc.spatial[, 4]),
                     fill = "#6baed680")

p <- p + geom_line(aes(x = 1:feature.count, y = auc.random[, 3]), color = "#99000d")
p <- p + geom_point(aes(x = 1:feature.count, y = auc.random[, 3]), color = "#99000d")

p <- p + geom_line(aes(x = 1:feature.count, y = auc.spatial[, 3]), color = "#084594")
p <- p + geom_point(aes(x = 1:feature.count, y = auc.spatial[, 3]), color = "#084594")

p <- p + geom_vline(xintercept = feature.random, linetype = "dotted")
p <- p + geom_point(aes(x = feature.random, y = auc.random[feature.random, 3]), size = 4)

p <- p + geom_vline(xintercept = feature.spatial, linetype = "dotted")
p <- p + geom_point(aes(x = feature.spatial, y = auc.spatial[feature.spatial, 3]), size = 4)

p <- p + scale_x_continuous(name = "# of Features",
                            limits = c(0, feature.count + 1),
                            breaks = seq(0, 100, by = 4),
                            expand = c(0, 0))
p <- p + scale_y_continuous(name = "AUC",
                            limits = c(min(auc.spatial[, 1]), 1.0),
#                           breaks = seq(0.0, 1.0, by = 0.02),
                            expand = c(0.01, 0.01))
p <- p + theme_bw()
p

ggsave(paste(s, "_AUC_Natural_O", n, "_", pool, "_", m, ".png", sep = ""),
       plot = p, path = s, width = 16, height = 10, units = "cm", dpi = 600)
