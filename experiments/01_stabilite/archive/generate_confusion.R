# ============================================================================
# Génère les matrices de confusion à partir des triplets stabilité-optimaux
# ============================================================================
# Usage : source("experiments/01_stabilite/generate_confusion.R")
# Prérequis : run_all_stability.R déjà exécuté (_best.csv présents)
# ============================================================================

library(cluster)
results_dir <- "experiments/01_stabilite/results"
DATASETS <- c("canadian", "aemet", "growth", "tecator")
preprocess_files <- c(
  canadian = "src/00_preprocess.R",
  aemet    = "src/00_preprocess_aemet.R",
  growth   = "src/00_preprocess_growth.R",
  tecator  = "src/00_preprocess_tecator.R"
)

set.seed(42)
for (d in DATASETS) {
  best <- read.csv(file.path(results_dir, paste0("stabilite_", d, "_best.csv")))
  source(preprocess_files[d])
  source("src/01_lissage.R")
  source("src/02_fpca.R")
  source("src/03_distances.R")

  a <- best$alpha
  w <- if ("omega" %in% names(best)) best$omega else best$r
  k_val <- best$k

  Dp <- sqrt((1 - a) * D0_norm^2 + a * D1_norm^2)
  Dw <- sqrt(w * Dp^2 + (1 - w) * Ds_norm^2)
  pam_best <- pam(as.dist(Dw), k = k_val, diss = TRUE)

  conf_mat <- as.matrix(table(regions, pam_best$clustering))
  write.csv(conf_mat, file.path(results_dir, paste0("confusion_", d, ".csv")))
}
cat("Matrices de confusion générées dans results/\n")
