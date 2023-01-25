ecdf.matrix <- function(feature.vector, reference.frame, value.frame) {
  percentile.matrix <- matrix(nrow = nrow(value.frame), ncol = length(feature.vector))
  for (i in 1:length(feature.vector)) {
    percentile.matrix[, i] <- ecdf(reference.frame[, feature.vector[i]])(value.frame[, feature.vector[i]])
  }
  percentile.matrix
}

mess.matrix <- function(feature.vector, reference.frame, value.frame, percentile.matrix) {
  mess <- matrix(nrow = nrow(value.frame), ncol = length(feature.vector))
  for (i in 1:length(feature.vector)) {
    p <- percentile.matrix[, i]
    p.1 <- p == 0.0
    p.2 <- p > 0.0 & p <= 0.5
    p.3 <- p > 0.5 & p < 1.0
    p.4 <- p == 1.0

    r.min <- min(reference.frame[, feature.vector[i]])
    r.max <- max(reference.frame[, feature.vector[i]])
    r.range <- r.max - r.min
    v <- value.frame[, feature.vector[i]]
  
    m <- mess[, i]
    m[p.1] <- ((v[p.1] - r.min) / r.range) * 100
    m[p.2] <- p[p.2] * 200
    m[p.3] <- (1 - p[p.3]) * 200
    m[p.4] <- ((r.max - v[p.4]) / r.range) * 100
    mess[, i] <- m
  }
  mess
}

pearson.frame <- function(feature.vector, reference.frame) {
  pearson.matrix <- cor(reference.frame[, feature.vector], method = "pearson")
  feature.matrix <- combn(feature.vector, 2)

  pearson.vector <- c()
  for (i in 1:ncol(feature.matrix)) {
    pearson.vector <- c(pearson.vector, pearson.matrix[feature.matrix[1, i],
                                                       feature.matrix[2, i]])
  }
  cor.frame <- data.frame(Feature.A = feature.matrix[1, ],
                          Feature.B = feature.matrix[2, ],
                          R = pearson.vector, stringsAsFactors = FALSE)
  cor.frame[order(cor.frame$R, decreasing = TRUE), ]
}

messcor <- function(feature.vector, mess.values, cor.frame) {
  feature.a <- match(cor.frame$Feature.A, feature.vector)
  feature.b <- match(cor.frame$Feature.B, feature.vector)

  mess.a.vector <- c()
  mess.b.vector <- c()
  remove.vector <- c()
  for (i in 1:nrow(cor.frame)) {
    mess.a <- mean(mess.values[, feature.a[i]])
    mess.b <- mean(mess.values[, feature.b[i]])
    r <- feature.vector[c(feature.a[i], feature.b[i])[which.min(c(mess.a, mess.b))]]

    mess.a.vector <- c(mess.a.vector, mess.a)
    mess.b.vector <- c(mess.b.vector, mess.b)
    remove.vector <- c(remove.vector, r)
  }
  messcor.frame <- data.frame(Feature.A = cor.frame$Feature.A,
                              Feature.B = cor.frame$Feature.B,
                              R = cor.frame$R,
                              MESS.A = mess.a.vector,
                              MESS.B = mess.b.vector,
                              Rem = remove.vector, stringsAsFactors = FALSE)
}

r.vector <- function(mess.frame, threshold) {
  r.feature.vector <- c()

  while (mess.frame[1, "R"] > threshold) {
    r.feature <- mess.frame[1, "Rem"]
    r.rows <- c()

    for (i in 1:nrow(mess.frame)) {
      if (mess.frame[i, "Feature.A"] == r.feature | mess.frame[i, "Feature.B"] == r.feature) {
        r.rows <- c(r.rows, i)
      }
    }

    r.feature.vector <- c(r.feature.vector, r.feature)
    mess.frame <- mess.frame[-1 * r.rows, ]
  }
  r.feature.vector
}

mess.stats <- function(set.name, feature.vector, selected.vector, mess.values, realms.unique, realm.vector) {
  stat.matrix <- matrix(nrow = 2, ncol = length(realms.unique) + 1)
  colnames(stat.matrix) <- c("Global", realms.unique)
  row.names(stat.matrix) <- c(paste(set.name, ".MESS", sep = ""), paste(set.name, ".Binary", sep = ""))
  
  values <- apply(mess.values[, match(selected.vector, feature.vector)], 1, min)
  binary <- values > 0
  
  stat.matrix[1, 1] <- mean(values)
  stat.matrix[2, 1] <- sum(binary) / length(binary)
  for (i in 1:length(realms.unique)) {
    stat.matrix[1, i + 1] <- mean(values[realm.vector == realms.unique[i]])
    stat.matrix[2, i + 1] <- sum(binary[realm.vector == realms.unique[i]]) / sum(realm.vector == realms.unique[i])
  }
  stat.matrix
}
