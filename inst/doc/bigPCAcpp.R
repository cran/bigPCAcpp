## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----bigmatrix-setup----------------------------------------------------------
library(bigmemory)
library(bigPCAcpp)

iris_mat <- as.matrix(iris[, 1:4])
big_iris <- as.big.matrix(iris_mat, type = "double")

## ----pca-run------------------------------------------------------------------
big_pca <- pca_bigmatrix(
  xpMat = big_iris,
  center = TRUE,
  scale = TRUE,
  ncomp = 4L,
  block_size = 128L
)
str(big_pca)

## ----prcomp-compare-----------------------------------------------------------
base_pca <- prcomp(iris_mat, center = TRUE, scale. = TRUE)

align_columns <- function(reference, target) {
  aligned <- target
  cols <- min(ncol(reference), ncol(target))
  for (j in seq_len(cols)) {
    ref <- reference[, j]
    tgt <- target[, j]
    if (sum((ref - tgt)^2) > sum((ref + tgt)^2)) {
      aligned[, j] <- -tgt
    }
  }
  aligned
}

rotation_aligned <- align_columns(base_pca$rotation, big_pca$rotation)
max_rotation_error <- max(abs(rotation_aligned - base_pca$rotation))
max_sdev_error <- max(abs(big_pca$sdev - base_pca$sdev))

big_scores <- pca_scores_bigmatrix(
  xpMat = big_iris,
  rotation = big_pca$rotation,
  center = big_pca$center,
  scale = big_pca$scale,
  block_size = 128L
)

scores_aligned <- align_columns(base_pca$x, big_scores)
max_score_error <- max(abs(scores_aligned - base_pca$x))

c(
  rotation = max_rotation_error,
  sdev = max_sdev_error,
  scores = max_score_error
)

## ----diagnostics--------------------------------------------------------------
loadings <- pca_variable_loadings(big_pca$rotation, big_pca$sdev)
correlations <- pca_variable_correlations(
  big_pca$rotation,
  big_pca$sdev,
  big_pca$column_sd,
  big_pca$scale
)
contributions <- pca_variable_contributions(loadings)

head(loadings)
head(correlations)
head(contributions)
range(correlations)

## ----plot-scree, fig.cap="Scree plot of variance explained by each component."----
pca_plot_scree(big_pca)

## ----plot-scores, fig.cap="Scores for the first two principal components."----
pca_plot_scores(
  big_iris,
  rotation = big_pca$rotation,
  center = big_pca$center,
  scale = big_pca$scale,
  max_points = nrow(big_iris),
  sample = "head"
)

## ----plot-correlation, fig.cap="Correlation circle highlighting how variables align with the first two components."----
pca_plot_correlation_circle(
  correlations,
  components = c(1L, 2L)
)

## ----plot-biplot, fig.cap="Biplot combining sample scores and variable loadings."----
pca_plot_biplot(
  big_scores,
  loadings,
  components = c(1L, 2L)
)

## ----svd-example--------------------------------------------------------------
svd_res <- svd_bigmatrix(big_iris, nu = 2L, nv = 2L, block_size = 128L)
svd_res$d

## ----robust-pca---------------------------------------------------------------
robust_pca <- pca_robust(iris_mat, ncomp = 4L)
robust_pca$sdev
robust_pca$robust_weights[1:10]

## ----robust-svd---------------------------------------------------------------
robust_svd <- svd_robust(iris_mat, ncomp = 3L)
robust_svd$d
robust_svd$weights[1:10]

## ----filebacked-example-------------------------------------------------------
library(bigmemory)
library(bigPCAcpp)

path <- tempfile(fileext = ".bin")
desc <- paste0(path, ".desc")

bm <- filebacked.big.matrix(
  nrow = nrow(iris_mat),
  ncol = ncol(iris_mat),
  type = "double",
  backingfile = basename(path),
  backingpath = dirname(path),
  descriptorfile = basename(desc)
)

bm[,] <- iris_mat

pca <- pca_bigmatrix(bm, center = TRUE, scale = TRUE, ncomp = 4)
scores <- filebacked.big.matrix(
  nrow = nrow(bm),
  ncol = ncol(pca$rotation),
  type = "double",
  backingfile = "scores.bin",
  backingpath = dirname(path),
  descriptorfile = "scores.desc"
)

pca_scores_stream_bigmatrix(
  bm,
  scores,
  pca$rotation,
  center = pca$center,
  scale = pca$scale
)

## ----filebacked-plot, fig.cap="Scores streamed from a file-backed big.matrix."----
pca_plot_scores(
  bm,
  rotation = pca$rotation,
  center = pca$center,
  scale = pca$scale,
  components = c(1L, 2L),
  max_points = nrow(bm),
  sample = "head"
)

## ----filebacked-example-2, eval = FALSE---------------------------------------
# library(bigmemory)
# library(bigPCAcpp)
# 
# path <- tempfile(fileext = ".bin")
# desc <- paste0(path, ".desc")
# 
# bm <- filebacked.big.matrix(
#   nrow = 5000,
#   ncol = 50,
#   type = "double",
#   backingfile = basename(path),
#   backingpath = dirname(path),
#   descriptorfile = basename(desc)
# )
# 
# pca <- pca_bigmatrix(bm, center = TRUE, scale = TRUE, ncomp = 5)
# scores <- filebacked.big.matrix(
#   nrow = nrow(bm),
#   ncol = ncol(pca$rotation),
#   type = "double",
#   backingfile = "scores.bin",
#   backingpath = dirname(path),
#   descriptorfile = "scores.desc"
# )
# 
# pca_scores_stream_bigmatrix(
#   bm,
#   scores,
#   pca$rotation,
#   center = pca$center,
#   scale = pca$scale
# )

