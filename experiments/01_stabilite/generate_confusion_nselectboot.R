# ============================================================================
# Matrices de confusion à partir des résultats nselectboot (Fang & Wang)
# ============================================================================
# Règle : premier couple (α, ω) où k_opt = k_vrai ; sinon k_opt le plus proche.
#
# Usage : source("experiments/01_stabilite/generate_confusion_nselectboot.R")
# Prérequis : nselectboot_*.csv dans results/
# ============================================================================

library(cluster)
library(mclust)

results_dir <- "experiments/01_stabilite/results"
DATASETS    <- c("canadian", "aemet", "growth", "tecator")
k_vrai      <- c(canadian = 4, aemet = 4, growth = 2, tecator = 3)
preprocess  <- c(
  canadian = "src/00_preprocess.R",
  aemet    = "src/00_preprocess_aemet.R",
  growth   = "src/00_preprocess_growth.R",
  tecator  = "src/00_preprocess_tecator.R"
)

set.seed(42)

cat("\n=== Matrices de confusion (nselectboot uniquement) ===\n\n")

for (d in DATASETS) {
  f_nb <- file.path(results_dir, paste0("nselectboot_", d, ".csv"))
  if (!file.exists(f_nb)) {
    cat(sprintf("  %s : nselectboot absent, skip\n", d))
    next
  }

  source(preprocess[d])
  source("src/01_lissage.R")
  source("src/02_fpca.R")
  source("src/03_distances.R")

  nb <- read.csv(f_nb)
  kv <- k_vrai[d]
  idx_ok <- which(nb$k_opt == kv)
  if (length(idx_ok) == 0) {
    idx_ok <- which.min(abs(nb$k_opt - kv))
  }
  idx_ok <- idx_ok[1]
  a_n <- nb$alpha[idx_ok]
  w_n <- if ("omega" %in% names(nb)) nb$omega[idx_ok] else nb$r[idx_ok]
  k_n <- nb$k_opt[idx_ok]

  Dp_n <- sqrt((1 - a_n) * D0_norm^2 + a_n * D1_norm^2)
  Dw_n <- sqrt(w_n * Dp_n^2 + (1 - w_n) * Ds_norm^2)
  pam_n <- pam(as.dist(Dw_n), k = k_n, diss = TRUE)

  conf_nb <- as.matrix(table(regions, pam_n$clustering))
  write.csv(conf_nb, file.path(results_dir, paste0("confusion_nselectboot_", d, ".csv")))
  ari_nb <- adjustedRandIndex(pam_n$clustering, as.integer(regions))

  cat(sprintf("  %s : nselectboot (k=%d, α=%.2f, ω=%.2f) ARI=%.3f\n",
    d, k_n, a_n, w_n, ari_nb))
}

cat("\n  Fichiers : confusion_nselectboot_*.csv\n\n")
