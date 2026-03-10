# ============================================================================
# Expérience 01 — Sélection par stabilité du clustering
# ============================================================================
#
# Protocole : pour chaque triplet (k, α, ω), on mesure la stabilité par
# bootstrap (sous-échantillonnage 80%, ARI entre partition complète et
# partition sur sous-échantillon). Le triplet avec stabilité max est retenu.
#
# Usage :
#   source("experiments/01_stabilite/run_stabilite.R")
# ou depuis R :
#   DATASET <- "canadian"
#   source("experiments/01_stabilite/run_stabilite.R")
#
# ============================================================================

library(cluster)
library(mclust)

# --- Paramètres (modifiables) ---
# Grille réduite pour prototype rapide ; augmenter pour résultats finals
ALPHAS <- c(0, 0.5, 1)      # prototype ; c(0, 0.1, ..., 1) pour grille fine
OMEGAS <- c(0, 0.5, 1)
K_VALUES <- 2:6
B <- 10                     # prototype ; 50 pour final
SUBSAMPLE_FRAC <- 0.8

# --- Dataset : utiliser DATASET si défini, sinon canadian ---
if (!exists("DATASET")) DATASET <- "canadian"
preprocess_file <- switch(DATASET,
  "canadian" = "src/00_preprocess.R",
  "tecator"  = "src/00_preprocess_tecator.R",
  "aemet"    = "src/00_preprocess_aemet.R",
  "growth"   = "src/00_preprocess_growth.R",
  "src/00_preprocess.R"
)

cat(sprintf("\n=== Stabilité du clustering — %s ===\n", DATASET))
cat(sprintf("  Grille : α ∈ {%s}, ω ∈ {%s}, k ∈ {%s}\n",
            paste(ALPHAS, collapse = ","), paste(OMEGAS, collapse = ","),
            paste(K_VALUES, collapse = ",")))
cat(sprintf("  B = %d, sous-échantillon = %.0f%%\n\n", B, SUBSAMPLE_FRAC * 100))

set.seed(42)

# --- Charger le pipeline jusqu'à l'étape 03 (distances) ---
source(preprocess_file)
source("src/01_lissage.R")
source("src/02_fpca.R")
source("src/03_distances.R")

n <- ncol(Y_brut)
n_subsample <- max(2, floor(n * SUBSAMPLE_FRAC))

# Vérité terrain (si disponible, pour évaluation finale)
labels_vrai <- if (exists("regions")) as.integer(regions) else NULL

# --- Grille des triplets ---
grille <- expand.grid(k = K_VALUES, alpha = ALPHAS, omega = OMEGAS)
grille$stabilite <- NA

# --- Boucle principale ---
cat("  Calcul de la stabilité...\n")
for (idx in 1:nrow(grille)) {
  k_val <- grille$k[idx]
  a <- grille$alpha[idx]
  w <- grille$omega[idx]

  # Distance Dw(α, ω)
  Dp <- sqrt((1 - a) * D0_norm^2 + a * D1_norm^2)
  Dw <- sqrt(w * Dp^2 + (1 - w) * Ds_norm^2)

  # Partition complète C
  pam_full <- pam(as.dist(Dw), k = k_val, diss = TRUE)
  C_full <- pam_full$clustering

  # Bootstrap
  aris <- numeric(B)
  for (b in 1:B) {
    idx_sub <- sort(sample(n, size = n_subsample, replace = FALSE))

    # Sous-matrice de distance
    Dw_sub <- Dw[idx_sub, idx_sub]

    # Partition sur sous-échantillon C_b*
    pam_sub <- pam(as.dist(Dw_sub), k = k_val, diss = TRUE)

    # C restreinte aux individus du sous-échantillon
    C_restricted <- C_full[idx_sub]

    # ARI entre les deux partitions (sur les mêmes individus)
    aris[b] <- adjustedRandIndex(C_restricted, pam_sub$clustering)
  }

  grille$stabilite[idx] <- mean(aris)
}

# --- Triplet optimal par stabilité ---
idx_best <- which.max(grille$stabilite)
k_best <- grille$k[idx_best]
alpha_best <- grille$alpha[idx_best]
omega_best <- grille$omega[idx_best]
stab_best <- grille$stabilite[idx_best]

cat(sprintf("\n  --- Résultat ---\n"))
cat(sprintf("  Triplet stabilité-optimal : k = %d, α = %.1f, ω = %.1f\n",
            k_best, alpha_best, omega_best))
cat(sprintf("  Stabilité = %.3f\n", stab_best))

# --- Si vérité terrain : ARI du triplet stabilité-optimal ---
if (!is.null(labels_vrai)) {
  Dp_best <- sqrt((1 - alpha_best) * D0_norm^2 + alpha_best * D1_norm^2)
  Dw_best <- sqrt(omega_best * Dp_best^2 + (1 - omega_best) * Ds_norm^2)
  pam_best <- pam(as.dist(Dw_best), k = k_best, diss = TRUE)
  ari_stab <- adjustedRandIndex(pam_best$clustering, labels_vrai)
  cat(sprintf("  ARI (vs vérité terrain, k=%d) = %.3f\n", k_best, ari_stab))
}

# --- Sauvegarde ---
results_dir <- "experiments/01_stabilite/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
fichier <- file.path(results_dir, sprintf("stabilite_%s.csv", DATASET))
write.csv(grille, fichier, row.names = FALSE)
cat(sprintf("\n  Résultats sauvegardés : %s\n", fichier))

# --- Top 5 triplets ---
grille_ord <- grille[order(-grille$stabilite), ]
cat("\n  Top 5 triplets (k, α, ω) par stabilité :\n")
print(head(grille_ord[, c("k", "alpha", "omega", "stabilite")], 5))
