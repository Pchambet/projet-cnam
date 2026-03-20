# ============================================================================
# ÉTAPE 00 (SIMULÉ) : PRÉPARATION DES DONNÉES HYBRIDES SIMULÉES
# ============================================================================
#
# Source : docs/biblio/notes/RE_Lectures_ACP_hybride/simulations.R (Cas2_deriv)
#
# Sorties alignées avec le pipeline existant :
# - Y_brut  : matrice (temps x individus) pour la partie fonctionnelle
# - Z       : data.frame (individus x variables) pour la partie vectorielle
# - regions : classes vraies (vérité terrain simulée)
# ============================================================================

library(fda)

cat("\n--- Étape 00 (simulé) : Préparation des données ---\n")

source("docs/biblio/notes/RE_Lectures_ACP_hybride/simulations.R")

# Reproductibilité
if (!exists("SEED_SIM")) SEED_SIM <- 42
set.seed(SEED_SIM)

# Paramètres simulés par défaut (modifiables avant source)
if (!exists("NC_SIM")) NC_SIM <- c(100, 100, 100)
if (!exists("LEN_T_SIM")) LEN_T_SIM <- 60
if (!exists("DELTA_SIM")) DELTA_SIM <- c(1, 1, 1)
if (!exists("SIGMA2_SIM")) SIGMA2_SIM <- 0.2
if (!exists("TAU2_SIM")) TAU2_SIM <- 0.2
# Dimension du bloc vectoriel Z (Cas2_deriv). Défaut 20 : ACP sur Z (02b) plus informative qu'avec p=10.
# Pour sensibilité / espace plus riche : assigner P_Z_SIM <- 30L avant source().
if (!exists("P_Z_SIM")) P_Z_SIM <- 20L
# Optionnel : nombre de composantes PCA sur Z après simulation (NULL = désactivé)
if (!exists("PCA_Z_Q")) PCA_Z_Q <- NULL

sim <- Cas2_deriv(
  nc = NC_SIM,
  len_t = LEN_T_SIM,
  delta = DELTA_SIM,
  p_z = as.integer(P_Z_SIM),
  sigma2 = SIGMA2_SIM,
  tau2 = TAU2_SIM
)

# Partie fonctionnelle observée:
# Cas2_deriv retourne X_obs en format (individus x temps).
# Le pipeline attend Y_brut en format (temps x individus).
Y_brut <- t(sim$X_obs)
t_jours <- 1:nrow(Y_brut)
n <- ncol(Y_brut)
N <- nrow(Y_brut)

# Partie vectorielle simulée
Z <- as.data.frame(sim$Y)
colnames(Z) <- paste0("z", seq_len(ncol(Z)))
noms_stations <- paste0("sim_", sprintf("%03d", seq_len(n)))
rownames(Z) <- noms_stations

# Optionnel : ACP classique sur Z (réduction de bruit / dimension avant 01–04)
if (!is.null(PCA_Z_Q) && is.numeric(PCA_Z_Q) && PCA_Z_Q > 0) {
  qz <- min(as.integer(PCA_Z_Q), ncol(Z))
  if (qz < ncol(Z)) {
    pz_obj <- prcomp(Z, center = TRUE, scale. = TRUE)
    Z <- as.data.frame(pz_obj$x[, seq_len(qz), drop = FALSE])
    colnames(Z) <- paste0("z", seq_len(qz))
    rownames(Z) <- noms_stations
    Z_pca_varprop <- summary(pz_obj)$importance[2, seq_len(qz)]
    cat(sprintf("  PCA(Z) : q = %d composantes (var. expl. cumulée ~ %.0f%%)\n",
                qz, 100 * sum(Z_pca_varprop)))
  } else {
    cat("  PCA(Z) : q >= p — ignoré.\n")
  }
}

# Vérité terrain
regions <- factor(paste0("C", sim$class), levels = paste0("C", sort(unique(sim$class))))

# Métadonnées pipeline
dataset_name <- "simulated_hybrid_cas2_deriv"
rangeval <- c(1, N)
xlab_courbe <- "Temps simulé"
ylab_courbe <- "Signal simulé"

cat(sprintf("  Dataset : %s\n", dataset_name))
cat(sprintf("  n = %d individus | N = %d points temporels\n", n, N))
cat(sprintf("  Variables Z : %d\n", ncol(Z)))
cat(sprintf("  Classes simulées : %s\n",
            paste(names(table(regions)), sprintf("(%d)", table(regions)), collapse = ", ")))
cat("  -> Prêt pour lissage (01), FPCA (02), ACP hybride (02b), noyaux hybrides (03b).\n")
