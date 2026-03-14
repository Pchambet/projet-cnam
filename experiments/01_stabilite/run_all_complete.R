# ============================================================================
# Exécution complète : stabilité + nselectboot + analyse + confusion + rapport
# ============================================================================
# Usage : Rscript experiments/01_stabilite/run_all_complete.R
# Ou : source("experiments/01_stabilite/run_all_complete.R")
# ============================================================================

# Se placer à la racine du projet
if (basename(getwd()) != "mars") {
  if (file.exists("experiments/01_stabilite/run_all_complete.R")) {
    # déjà à la racine
  } else if (file.exists("../experiments/01_stabilite/run_all_complete.R")) {
    setwd("..")
  } else {
    setwd("/Users/pierre/Desktop/mars")
  }
}

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  EXÉCUTION COMPLÈTE — Expérience 01\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# 1. Stabilité
cat(">>> Étape 1/5 : Stabilité (4 datasets, B=150, grille 21×21)\n")
source("experiments/01_stabilite/run_all_stability.R", local = FALSE)
cat("\n")

# 2. nselectboot
cat(">>> Étape 2/5 : nselectboot (4 datasets, B=150, grille 21×21)\n")
source("experiments/01_stabilite/run_all_nselectboot.R", local = FALSE)
cat("\n")

# 3. Heatmaps stabilité
cat(">>> Étape 3/5 : Heatmaps stabilité\n")
source("experiments/01_stabilite/analyse_stabilite.R", local = FALSE)
cat("\n")

# 4. Heatmaps nselectboot
cat(">>> Étape 4/5 : Heatmaps nselectboot\n")
source("experiments/01_stabilite/analyse_nselectboot.R", local = FALSE)
cat("\n")

# 5. Matrices de confusion (les deux méthodes)
cat(">>> Étape 5/5 : Matrices de confusion (stabilité + nselectboot)\n")
source("experiments/01_stabilite/generate_confusion_both.R", local = FALSE)
cat("\n")

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TERMINÉ. Compiler le rapport : make (depuis experiments/01_stabilite/)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
