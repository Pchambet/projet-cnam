# ============================================================================
# EXPÉRIENCE 03 : CHAÎNE ACP HYBRIDE + NOYAUX SUR DONNÉES SIMULÉES
# ============================================================================
#
# Pipeline exécuté :
#   00_preprocess_simulated.R -> 01_lissage.R -> 02_fpca.R
#   -> 02b_pca_hybride_reconstruction.R -> 03b_distances_noyaux_hybrides.R
# Nomenclature : r = ratio HFV (02b) ; ω = Dw (03) — voir STATE_OF_PROJECT.md
#
# Sortie principale :
#   ARI final sur PAM(k = nb classes simulées) avec D_K (noyaux hybrides)
# ============================================================================

source("setup.R")

set.seed(42)
SEED_SIM <- 42

# Paramètres simulation (ajustables avant source)
NC_SIM <- c(100, 100, 100)
LEN_T_SIM <- 60
DELTA_SIM <- c(1, 1, 1)
SIGMA2_SIM <- 0.2
TAU2_SIM <- 0.2

source("src/00_preprocess_simulated.R")
source("src/01_lissage.R")
source("src/02_fpca.R")
source("src/02b_pca_hybride_reconstruction.R")
source("src/03b_distances_noyaux_hybrides.R")

n_clusters <- length(levels(regions))
pam_sim <- cluster::pam(as.dist(D_K), k = n_clusters, diss = TRUE)
ari_sim <- mclust::adjustedRandIndex(pam_sim$clustering, as.integer(regions))
sil_sim <- summary(cluster::silhouette(pam_sim$clustering, as.dist(D_K)))$avg.width

cat("\n=== RESULTATS SIMULES (ACP hybride + noyaux hybrides) ===\n")
cat(sprintf("dim(D_K) = %d x %d\n", nrow(D_K), ncol(D_K)))
cat(sprintf("NA dans D_K = %d\n", sum(is.na(D_K))))
cat(sprintf("k = %d\n", n_clusters))
cat(sprintf("Silhouette = %.6f\n", sil_sim))
cat(sprintf("ARI final = %.6f\n", ari_sim))

cat("\n=== MATRICE DE CONFUSION (clusters vs classes simulées) ===\n")
print(table(cluster = pam_sim$clustering, classe = regions))
