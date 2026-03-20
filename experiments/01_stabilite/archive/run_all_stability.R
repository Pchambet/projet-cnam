# ============================================================================
# Expérience 01 — Exécution stabilité sur les 4 datasets
# ============================================================================
#
# OBSOLÈTE — déplacé dans archive/. Usage (depuis la racine du projet) :
#   source("experiments/01_stabilite/archive/run_all_stability.R")
#
# ============================================================================

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  EXPÉRIENCE 01 — SÉLECTION PAR STABILITÉ (4 datasets)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# --- Étape 1 : Paramétrage global ---
DATASETS <- c("canadian", "aemet", "growth", "tecator")

# Grille fine, B élevé pour robustesse
ALPHAS   <- seq(0, 1, by = 0.05)   # 21×21 grille (pas 0.05)
OMEGAS   <- seq(0, 1, by = 0.05)   # grain fin sur [0,1]
K_VALUES <- 2:6                    # k = 2, 3, 4, 5, 6
B        <- 150                    # répétitions bootstrap
# Sous-échantillon : 80 % = compromis standard. Plus haut (90 %) = moins de
# perturbation, ARI plus stables mais moins discriminatif. Plus bas (63 %, 50 %)
# = plus de variabilité, mieux pour détecter l'instabilité mais risque de
# sur-pénaliser les grands k quand n est petit. Avec n=35 (Canadian), 80 % →
# 28 individus, suffisant pour PAM à k≤6.
SUBSAMPLE_FRAC <- 0.8              # 80 % des individus par réplication

cat("Paramètres :\n")
cat(sprintf("  Datasets   : %s\n", paste(DATASETS, collapse = ", ")))
cat(sprintf("  α          : %s\n", paste(round(ALPHAS, 2), collapse = ", ")))
cat(sprintf("  ω          : %s\n", paste(round(OMEGAS, 2), collapse = ", ")))
cat(sprintf("  k          : %s\n", paste(K_VALUES, collapse = ", ")))
cat(sprintf("  B          : %d\n", B))
cat(sprintf("  Sous-échantillon : %.0f%%\n\n", SUBSAMPLE_FRAC * 100))

set.seed(42)

# --- Étape 2 : Exécution pour chaque dataset ---
for (d in DATASETS) {
  DATASET <- d
  source("experiments/01_stabilite/archive/run_stabilite.R", local = FALSE)
}

cat("\n══════════════════════════════════════════════════════════════════════\n")
cat("  Exécution terminée. Résultats dans experiments/01_stabilite/results/\n")
cat("  Lancer : source(\"experiments/01_stabilite/archive/analyse_stabilite.R\")\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
