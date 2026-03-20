# ============================================================================
# Expérience 01 — Sélection de k par nselectboot (Fang & Wang 2012)
# ============================================================================
#
# nselectboot minimise l'instabilité (pas maximise la stabilité ARI).
# Méthode : paires de bootstrap, comparaison des partitions après
# classification croisée. k optimal = argmin Inst(k).
#
# Usage :
#   DATASET <- "canadian"
#   source("experiments/01_stabilite/run_nselectboot.R")
#
# ============================================================================

library(cluster)
library(mclust)
library(fpc)

if (!exists("DATASET")) DATASET <- "canadian"
# Valeurs par défaut : B=150, grille 21×21
if (!exists("B_NSELECT")) B_NSELECT <- 150
if (!exists("ALPHAS_NB")) ALPHAS_NB <- seq(0, 1, by = 0.05)
if (!exists("OMEGAS_NB")) OMEGAS_NB <- seq(0, 1, by = 0.05)
if (!exists("K_RANGE")) K_RANGE <- 2:6

preprocess_file <- switch(DATASET,
  "canadian" = "src/00_preprocess.R",
  "tecator"  = "src/00_preprocess_tecator.R",
  "aemet"    = "src/00_preprocess_aemet.R",
  "growth"   = "src/00_preprocess_growth.R",
  "src/00_preprocess.R"
)

cat(sprintf("\n=== nselectboot (Fang-Wang) — %s ===\n", DATASET))
cat(sprintf("  B = %d, k ∈ {%s}\n", B_NSELECT, paste(K_RANGE, collapse = ",")))
cat(sprintf("  Grille (α, ω) : %d × %d points\n\n", length(ALPHAS_NB), length(OMEGAS_NB)))

set.seed(42)

# Charger données et distances
source(preprocess_file)
source("src/01_lissage.R")
source("src/02_fpca.R")
source("src/03_distances.R")

n <- ncol(Y_brut)
labels_vrai <- if (exists("regions")) as.integer(regions) else NULL
k_vrai <- if (!is.null(labels_vrai)) length(unique(labels_vrai)) else NA

results_dir <- "experiments/01_stabilite/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# Grille (α, ω) pour tester
grille_nb <- expand.grid(alpha = ALPHAS_NB, omega = OMEGAS_NB)
grille_nb$k_opt <- NA
grille_nb$instability_min <- NA
grille_nb$stabk <- NA  # sera une liste ou JSON des stabk par k

all_stabk <- list()

for (idx in 1:nrow(grille_nb)) {
  a <- grille_nb$alpha[idx]
  w <- grille_nb$omega[idx]

  Dp <- sqrt((1 - a) * D0_norm^2 + a * D1_norm^2)
  Dw <- sqrt(w * Dp^2 + (1 - w) * Ds_norm^2)
  D_dist <- as.dist(Dw)

  cat(sprintf("  (α=%.1f, ω=%.1f) ", a, w))

  # nselectboot : clustermethod=claraCBI (PAM), classification="centroid"
  # Pour dist, claraCBI avec usepam=TRUE utilise PAM. Les médoïdes servent de centroïdes.
  res <- tryCatch({
    nselectboot(D_dist,
      B = B_NSELECT,
      distances = TRUE,
      clustermethod = claraCBI,
      usepam = TRUE,
      classification = "centroid",
      centroidname = NULL,  # reconnu auto pour claraCBI
      krange = K_RANGE,
      count = FALSE
    )
  }, error = function(e) {
    warning(sprintf("nselectboot erreur (α=%.1f, ω=%.1f): %s", a, w, e$message))
    NULL
  })

  if (is.null(res)) {
    cat("ERREUR\n")
    next
  }

  grille_nb$k_opt[idx] <- res$kopt
  grille_nb$instability_min[idx] <- min(res$stabk)
  grille_nb$stabk[idx] <- paste(round(res$stabk, 4), collapse = ";")
  cat(sprintf("k_opt = %d (instab. min = %.4f)\n", res$kopt, min(res$stabk)))
}

# Résumé
cat("\n--- Résumé nselectboot ---\n")
print(grille_nb[, c("alpha", "omega", "k_opt", "instability_min")])

# Comparaison avec vérité terrain
if (!is.na(k_vrai)) {
  cat(sprintf("\n  k_vrai = %d\n", k_vrai))
  cat(sprintf("  k_nselectboot (médiane sur grille) = %d\n",
    as.integer(median(grille_nb$k_opt, na.rm = TRUE))))
}

# Sauvegarde
fichier <- file.path(results_dir, sprintf("nselectboot_%s.csv", DATASET))
write.csv(grille_nb, fichier, row.names = FALSE)
cat(sprintf("\n  Résultats sauvegardés : %s\n", fichier))
