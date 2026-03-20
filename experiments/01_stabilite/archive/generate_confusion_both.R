# ============================================================================
# Génère les matrices de confusion pour les DEUX méthodes
# ============================================================================
# Stabilité : triplet (k, α, ω) du meilleur selon stabilité
# nselectboot : triplet où k_opt = k_vrai (ou le plus proche)
#
# Usage : source("experiments/01_stabilite/generate_confusion_both.R")
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

cat("\n=== Matrices de confusion (stabilité + nselectboot) ===\n\n")

for (d in DATASETS) {
  source(preprocess[d])
  source("src/01_lissage.R")
  source("src/02_fpca.R")
  source("src/03_distances.R")

  # --- Stabilité ---
  best_stab <- read.csv(file.path(results_dir, paste0("stabilite_", d, "_best.csv")))
  a_s <- best_stab$alpha
  w_s <- if ("omega" %in% names(best_stab)) best_stab$omega else best_stab$r
  k_s <- best_stab$k

  Dp_s <- sqrt((1 - a_s) * D0_norm^2 + a_s * D1_norm^2)
  Dw_s <- sqrt(w_s * Dp_s^2 + (1 - w_s) * Ds_norm^2)
  pam_s <- pam(as.dist(Dw_s), k = k_s, diss = TRUE)

  conf_stab <- as.matrix(table(regions, pam_s$clustering))
  write.csv(conf_stab, file.path(results_dir, paste0("confusion_stabilite_", d, ".csv")))
  ari_stab <- adjustedRandIndex(pam_s$clustering, as.integer(regions))

  # --- nselectboot : prendre (α, ω) où k_opt = k_vrai, ou le plus proche
  f_nb <- file.path(results_dir, paste0("nselectboot_", d, ".csv"))
  if (!file.exists(f_nb)) {
    cat(sprintf("  %s : nselectboot absent, skip\n", d))
    next
  }

  nb <- read.csv(f_nb)
  kv <- k_vrai[d]
  # Chercher une ligne où k_opt = k_vrai
  idx_ok <- which(nb$k_opt == kv)
  if (length(idx_ok) == 0) {
    # Sinon prendre celle où k_opt est le plus proche de k_vrai
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

  cat(sprintf("  %s : stab (k=%d, α=%.2f, ω=%.2f) ARI=%.3f | nselectboot (k=%d, α=%.2f, ω=%.2f) ARI=%.3f\n",
    d, k_s, a_s, w_s, ari_stab, k_n, a_n, w_n, ari_nb))
}

cat("\n  Fichiers : confusion_stabilite_*.csv, confusion_nselectboot_*.csv\n\n")
