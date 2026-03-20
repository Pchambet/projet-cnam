# ============================================================================
# Exécution complète — Expérience 01 (instabilité / nselectboot uniquement)
# ============================================================================
# nselectboot → heatmaps instabilité → matrices de confusion
#
# Usage : Rscript experiments/01_stabilite/run_all_complete.R
# Ou : source("experiments/01_stabilite/run_all_complete.R")
# ============================================================================

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
cat("  EXÉCUTION COMPLÈTE — Expérience 01 (instabilité / nselectboot)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat(">>> Étape 1/3 : nselectboot (4 datasets, B=150, grille 21×21)\n")
source("experiments/01_stabilite/run_all_nselectboot.R", local = FALSE)
cat("\n")

cat(">>> Étape 2/3 : Heatmaps instabilité (nselectboot)\n")
source("experiments/01_stabilite/analyse_nselectboot.R", local = FALSE)
cat("\n")

cat(">>> Étape 3/3 : Matrices de confusion (nselectboot)\n")
source("experiments/01_stabilite/generate_confusion_nselectboot.R", local = FALSE)
cat("\n")

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TERMINÉ. Compiler le rapport : make (depuis experiments/01_stabilite/)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
