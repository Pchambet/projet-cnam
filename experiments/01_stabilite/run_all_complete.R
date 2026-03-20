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

# Exécution complète : inclure le volet simulé après les 4 jeux réels (voir protocole.md)
if (!exists("RUN_NSELECTBOOT_SIMULATED")) RUN_NSELECTBOOT_SIMULATED <- TRUE

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  EXÉCUTION COMPLÈTE — Expérience 01 (instabilité / nselectboot)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

cat(">>> Étape 1/4 : nselectboot (4 datasets + volet simulé, B=150, grille 21×21)\n")
source("experiments/01_stabilite/run_all_nselectboot.R", local = FALSE)
cat("\n")

cat(">>> Étape 2/4 : Heatmaps instabilité (nselectboot, 4 jeux réels)\n")
source("experiments/01_stabilite/analyse_nselectboot.R", local = FALSE)
cat("\n")

cat(">>> Étape 3/4 : Heatmaps instabilité (nselectboot, données simulées)\n")
source("experiments/01_stabilite/analyse_nselectboot_simulated.R", local = FALSE)
cat("\n")

cat(">>> Étape 4/4 : Matrices de confusion (nselectboot, jeux réels)\n")
source("experiments/01_stabilite/generate_confusion_nselectboot.R", local = FALSE)
cat("\n")

cat("══════════════════════════════════════════════════════════════════════\n")
cat("  TERMINÉ. Compiler le rapport : make (depuis experiments/01_stabilite/)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
