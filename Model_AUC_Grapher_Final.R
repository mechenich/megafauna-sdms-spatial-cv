library(ggplot2)
z <- "/home/ad/home/m/mechenic/Projects/megafauna-sdms-spatial-cv"
setwd(z)

neighbors <- c("30", "20", "30", "06", "15", "20", "20", "11", "06", "04")

specieslist <- c("Alces_Alces", "Bison_Bonasus", "Camelus_Dromedarius",
                 "Ceratotherium_Simum", "Cervus_Elaphus", "Diceros_Bicornis",
                 "Elephas_Maximus", "Equus_Africanus", "Ovibos_Moschatus",
                 "Rhinoceros_Unicornis")

auc <- data.frame(species = character(),
                  experiment = character(),
                  pool = character(),
                  validation = character(),
                  auc.0 = double(), auc.25 = double(), auc.50 = double(),
                  auc.75 = double(), auc.100 = double())

for (s in specieslist) {
    for (p in c("VIF9", "MESS9")) {
        for (cv in c("R10", "S10")) {
            load(paste(s, "/", s, "_Features_Natural_O",
                       neighbors[match(s, specieslist)], "_", cv, "_", p,
                       "_GLM.rdata", sep = ""))

            ffs.directory <- paste(s, "_Natural_O", neighbors[match(s, specieslist)],
                                   "_", cv, "_", p, "_GLM_AUC/", sep = "")
            ffs.files <- list.files(ffs.directory)

            e.auc <- double()            
            for (ffs in ffs.files) {
                load(paste(ffs.directory, ffs, sep = ""))
                e.auc <- c(e.auc, results$AUC.Valid[results$Feature.Index == nrow(f.set)])
            }

            ffs.quart <- quantile(e.auc, probs = seq(0, 1, 0.25))
            
            # -----------------------------------------------------------------
            load(paste(s, "/", s, "_AUC_Natural_O",
                       neighbors[match(s, specieslist)], "_", cv, "_", p,
                       "_GLM.rdata", sep = ""))

            fin.quart <- quantile(e.auc, probs = seq(0, 1, 0.25))

            e.frame <- data.frame(species = c(s, s),
                                  experiment = c("FFS", "Final"),
                                  pool = c(p, p),
                                  validation = c(cv, cv),
                                  auc.0 = c(ffs.quart[1], fin.quart[1]),
                                  auc.25 = c(ffs.quart[2], fin.quart[2]),
                                  auc.50 = c(ffs.quart[3], fin.quart[3]),
                                  auc.75 = c(ffs.quart[4], fin.quart[4]),
                                  auc.100 = c(ffs.quart[5], fin.quart[5]))
            auc <- rbind(auc, e.frame)
        }
    }
}

# -----------------------------------------------------------------------------
ffs.fin <- data.frame(x = auc$auc.50[auc$experiment == "FFS" & auc$validation == "S10"],
                      xmin = auc$auc.25[auc$experiment == "FFS" & auc$validation == "S10"],
                      xmax = auc$auc.75[auc$experiment == "FFS" & auc$validation == "S10"],
                      y = auc$auc.50[auc$experiment == "Final" & auc$validation == "S10"],
                      ymin = auc$auc.25[auc$experiment == "Final" & auc$validation == "S10"],
                      ymax = auc$auc.75[auc$experiment == "Final" & auc$validation == "S10"])
summary(lm(y ~ x, ffs.fin))

p <- ggplot()
p <- p + geom_errorbar(aes(x = ffs.fin$x, ymin = ffs.fin$ymin, ymax = ffs.fin$ymax),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")
p <- p + geom_errorbar(aes(xmin = ffs.fin$xmin, xmax = ffs.fin$xmax, y = ffs.fin$y),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")

p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#4c4c4c")
p <- p + geom_abline(intercept = -0.04314, slope = 1.04134, color = "#4c4c4c")

p <- p + scale_x_continuous("FFS AUC", limits = c(0.7, 1), expand = c(0, 0))
p <- p + scale_y_continuous("Final AUC", limits = c(0.7, 1), expand = c(0, 0))

p <- p + geom_point(aes(x = ffs.fin$x, y = ffs.fin$y))
p <- p + theme_bw()
p

ggsave("AUC_GLM_FFS_Final.png", plot = p, path = z, width = 10, height = 10, units = "cm", dpi = 600)

# -----------------------------------------------------------------------------
r.s <- data.frame(x = auc$auc.50[auc$experiment == "Final" & auc$validation == "R10"],
                  xmin = auc$auc.25[auc$experiment == "Final" & auc$validation == "R10"],
                  xmax = auc$auc.75[auc$experiment == "Final" & auc$validation == "R10"],
                  y = auc$auc.50[auc$experiment == "Final" & auc$validation == "S10"],
                  ymin = auc$auc.25[auc$experiment == "Final" & auc$validation == "S10"],
                  ymax = auc$auc.75[auc$experiment == "Final" & auc$validation == "S10"])
summary(lm(y ~ x, r.s))

p <- ggplot()
p <- p + geom_errorbar(aes(x = r.s$x, ymin = r.s$ymin, ymax = r.s$ymax),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")
p <- p + geom_errorbar(aes(xmin = r.s$xmin, xmax = r.s$xmax, y = r.s$y),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")

p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#4c4c4c")
p <- p + geom_abline(intercept = 0.34516, slope = 0.63866, color = "#4c4c4c")

p <- p + scale_x_continuous("R-CV Sets AUC", limits = c(0.7, 1), expand = c(0, 0))
p <- p + scale_y_continuous("S-CV Sets AUC", limits = c(0.7, 1), expand = c(0, 0))

p <- p + geom_point(aes(x = r.s$x, y = r.s$y))
p <- p + theme_bw()
p

ggsave("AUC_GLM_R10_S10.png", plot = p, path = z, width = 10, height = 10, units = "cm", dpi = 600)

# -----------------------------------------------------------------------------
vif.mess <- data.frame(x = auc$auc.50[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "VIF9"],
                       xmin = auc$auc.25[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "VIF9"],
                       xmax = auc$auc.75[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "VIF9"],
                       y = auc$auc.50[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "MESS9"],
                       ymin = auc$auc.25[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "MESS9"],
                       ymax = auc$auc.75[auc$experiment == "Final" & auc$validation == "S10" & auc$pool == "MESS9"])
summary(lm(y ~ x, vif.mess))

p <- ggplot()
p <- p + geom_errorbar(aes(x = vif.mess$x, ymin = vif.mess$ymin, ymax = vif.mess$ymax),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")
p <- p + geom_errorbar(aes(xmin = vif.mess$xmin, xmax = vif.mess$xmax, y = vif.mess$y),
                       width = 0.0045, linewidth = 0.45, color = "#bfbfbf")

p <- p + geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#4c4c4c")
p <- p + geom_abline(intercept = 0.01033, slope = 0.98864, color = "#4c4c4c")

p <- p + scale_x_continuous("VIF Sets AUC", limits = c(0.7, 1), expand = c(0, 0))
p <- p + scale_y_continuous("MESS Sets AUC", limits = c(0.7, 1), expand = c(0, 0))

p <- p + geom_point(aes(x = vif.mess$x, y = vif.mess$y))
p <- p + theme_bw()
p

ggsave("AUC_GLM_VIF9_MESS9.png", plot = p, path = z, width = 10, height = 10, units = "cm", dpi = 600)
