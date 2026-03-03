# ============================================================================
# Pipeline principal — Classification non-supervisée de données mixtes
# Stage : Application en océanographie
# Auteur : Pierre
# Date : Mars 2026
# ============================================================================
#
# Ce script exécute l'ensemble du pipeline :
#   1. Pré-traitement des données
#   2. Lissage (B-splines pénalisées)
#   3. ACP fonctionnelle
#   4. Distances (D0, D1, Dω, DK)
#   5. Clustering (PAM) et évaluation
#   6. Visualisation
#
# Usage : source("src/main.R")
# ============================================================================

cat("=== Pipeline FDA — Données mixtes ===\n\n")

# --- Reproductibilité ---
set.seed(42)

# --- Chargement des scripts ---
# source("src/00_preprocess.R")
# source("src/01_lissage.R")
# source("src/02_fpca.R")
# source("src/03_distances.R")
# source("src/04_clustering.R")
# source("src/05_visualisation.R")

cat("\n=== Pipeline terminé ===\n")
