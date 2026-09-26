# Synthetic CosMx-shaped data so the explorer runs without the real object.
# 126k cells, 4 samples, 2 groups, 3 clustering resolutions, auto + manual annotation.
set.seed(1)

N   <- 126000
K   <- 16                                   # "true" subpopulations
SAM <- c("HC_1", "HC_2", "SLE_1", "SLE_2")
GRP <- c(HC_1 = "HC", HC_2 = "HC", SLE_1 = "SLE", SLE_2 = "SLE")

sample_id <- factor(sample(SAM, N, TRUE), levels = SAM)
k         <- sample(K, N, TRUE)

# UMAP: one gaussian blob per subpopulation
ck <- cbind(runif(K, -11, 11), runif(K, -11, 11))
u1 <- ck[k, 1] + rnorm(N, 0, .9)
u2 <- ck[k, 2] + rnorm(N, 0, .9)

# Spatial: each subpopulation sits in its own domain, per sample
dk <- cbind(runif(K, 400, 4600), runif(K, 400, 4600))
x  <- dk[k, 1] + rnorm(N, 0, 260)
y  <- dk[k, 2] + rnorm(N, 0, 260)

# Resolutions are SHUFFLED, not nested-in-order: cluster 3 at res 0.3 has nothing
# to do with cluster 3 at res 0.8. That is how real re-clustering behaves.
relabel <- function(v) { m <- sample(sort(unique(v))); factor(m[match(v, sort(unique(v)))]) }
res_0.8 <- relabel(k)
res_0.5 <- relabel(ceiling(k / 2))
res_0.3 <- relabel(ceiling(k / 4))

types <- c("Podocyte", "PT", "TAL", "DCT", "Endothelial",
           "Fibroblast", "T cell", "B cell")
anno_auto <- factor(types[((k - 1) %% length(types)) + 1], levels = types)

# manual differs from auto: PT split by segment, plus a relabelled slice
anno_manual <- as.character(anno_auto)
anno_manual[anno_auto == "PT" & k > 8] <- "PT_S3"
anno_manual[anno_auto == "T cell" & runif(N) < .18] <- "NK"
anno_manual <- factor(anno_manual)

cells <- data.frame(
  cell_id = sprintf("c%06d", seq_len(N)),
  x, y, UMAP_1 = u1, UMAP_2 = u2,
  PCA_1 = u1 * .6 + rnorm(N, 0, 2), PCA_2 = u2 * .6 + rnorm(N, 0, 2),
  sample = sample_id, group = factor(GRP[as.character(sample_id)], levels = c("HC", "SLE")),
  res_0.3, res_0.5, res_0.8, anno_auto, anno_manual
)

# 12-vertex cell boundaries; i indexes back into cells (cheaper than a join)
V   <- 12
ang <- seq(0, 2 * pi, length.out = V + 1)[-(V + 1)]
rad <- 7 + runif(N, 0, 3)
i   <- rep(seq_len(N), each = V)
poly <- data.frame(i, px = x[i] + rad[i] * cos(ang), py = y[i] + rad[i] * sin(ang))

saveRDS(cells, "demo_cells.rds")
saveRDS(poly,  "demo_poly.rds",  compress = "xz")
cat(nrow(cells), "cells,", nrow(poly), "vertices\n")
