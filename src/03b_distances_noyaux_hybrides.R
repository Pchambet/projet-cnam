# ============================================================================
# ÉTAPE 03b : DISTANCES À NOYAUX HYBRIDES SUR COURBES RECONSTRUITES
# ============================================================================
#
# Entrées attendues :
#   - X_recon_mat : matrice (temps x individus) issue de 02b (reconstruction L²)
#   - Z           : tableau vectoriel **original** (même individus) — D_Z et K_s sur Z_std ;
#                   cohérent avec le reste du pipeline même si 02b fusionne des scores internes.
#
# ÉTAPE 1 : D_recon (L2) via intégration vectorisée
# ÉTAPE 2 : D_Z (euclidienne) sur Z standardisé
# ÉTAPE 3 : sigmas = médiane des distances non nulles
# ÉTAPE 4 : K_f, K_s, puis K_hybride = K_f * K_s
# ÉTAPE 5 : D_K depuis K_hybride, avec stabilisation numérique
# ============================================================================

library(cluster)
library(mclust)

cat("\n--- Étape 03b : Distances à noyaux hybrides ---\n")

.resolve_preprocess_file <- function(dataset) {
  switch(dataset,
    "canadian" = "src/00_preprocess.R",
    "tecator"  = "src/00_preprocess_tecator.R",
    "aemet"    = "src/00_preprocess_aemet.R",
    "growth"   = "src/00_preprocess_growth.R",
    stop(sprintf("Dataset inconnu '%s'.", dataset))
  )
}

.ensure_inputs_03b <- function() {
  if (!exists("DATASET")) DATASET <<- "canadian"

  if (!exists("Z") || !exists("regions") || !exists("t_jours")) {
    source(.resolve_preprocess_file(DATASET))
  }
  if (!exists("X_hat")) source("src/01_lissage.R")
  if (!exists("scores_fpca")) source("src/02_fpca.R")
  if (!exists("X_recon_mat")) source("src/02b_pca_hybride_reconstruction.R")
}

run_hybrid_kernel_distances <- function(verbose = TRUE) {
  .ensure_inputs_03b()

  if (ncol(X_recon_mat) != nrow(Z)) {
    stop(sprintf(
      "Incohérence dimensions : ncol(X_recon_mat)=%d vs nrow(Z)=%d",
      ncol(X_recon_mat), nrow(Z)
    ))
  }

  # --------------------------------------------------------------------------
  # ETAPE 1 : Distance L2 sur courbes reconstruites
  # --------------------------------------------------------------------------
  # D_recon_sq = ||x_i||² + ||x_j||² - 2<x_i,x_j>, avec intégration discrète.
  delta <- (max(t_jours) - min(t_jours)) / (length(t_jours) - 1)
  ss <- colSums(X_recon_mat^2) * delta
  inner <- crossprod(X_recon_mat) * delta
  D_recon_sq_raw <- outer(ss, ss, "+") - 2 * inner
  D_recon_sq <- pmax(0, zapsmall(D_recon_sq_raw))
  dim(D_recon_sq) <- dim(D_recon_sq_raw)
  dimnames(D_recon_sq) <- dimnames(D_recon_sq_raw)
  D_recon <- sqrt(D_recon_sq)

  # --------------------------------------------------------------------------
  # ETAPE 2 : Distance euclidienne sur Z standardisé
  # --------------------------------------------------------------------------
  Z_std <- scale(Z, center = TRUE, scale = TRUE)
  D_Z <- as.matrix(dist(Z_std))

  # --------------------------------------------------------------------------
  # ETAPE 3 : Echelles noyaux
  # --------------------------------------------------------------------------
  sigma_f <- median(D_recon[D_recon > 0])
  sigma_s <- median(D_Z[D_Z > 0])

  if (!is.finite(sigma_f) || sigma_f <= 0) {
    stop(sprintf("sigma_f invalide: %s", as.character(sigma_f)))
  }
  if (!is.finite(sigma_s) || sigma_s <= 0) {
    stop(sprintf("sigma_s invalide: %s", as.character(sigma_s)))
  }

  # --------------------------------------------------------------------------
  # ETAPE 4 : Noyau hybride
  # --------------------------------------------------------------------------
  K_f <- exp(-(D_recon^2) / (2 * sigma_f^2))
  K_s <- exp(-(D_Z^2) / (2 * sigma_s^2))
  K_hybride <- K_f * K_s

  # --------------------------------------------------------------------------
  # ETAPE 5 : Distance finale depuis le noyau hybride
  # --------------------------------------------------------------------------
  D_K_sq_raw <- outer(diag(K_hybride), diag(K_hybride), "+") - 2 * K_hybride
  D_K_sq <- pmax(0, zapsmall(D_K_sq_raw))
  dim(D_K_sq) <- dim(D_K_sq_raw)
  dimnames(D_K_sq) <- dimnames(D_K_sq_raw)
  D_K <- sqrt(D_K_sq)
  # Stabilisation finale numerique : symetrisation et diagonale exacte a 0.
  D_K <- as.matrix((D_K + t(D_K)) / 2)
  diag(D_K) <- 0

  # Contrôles qualité demandés
  sym_err <- max(abs(D_K - t(D_K)))
  diag_err <- max(abs(diag(D_K)))
  if (sym_err >= 1e-10) {
    stop(sprintf("D_K non symétrique (max err = %.3e)", sym_err))
  }
  if (diag_err >= 1e-10) {
    stop(sprintf("Diagonale D_K non nulle (max err = %.3e)", diag_err))
  }

  if (!is.null(rownames(Z))) {
    rownames(D_recon) <- colnames(D_recon) <- rownames(Z)
    rownames(D_Z) <- colnames(D_Z) <- rownames(Z)
    rownames(K_hybride) <- colnames(K_hybride) <- rownames(Z)
    rownames(D_K) <- colnames(D_K) <- rownames(Z)
  }

  if (verbose) {
    cat(sprintf("  sigma_f = %.6f | sigma_s = %.6f\n", sigma_f, sigma_s))
    cat(sprintf("  dim(D_K) = %d x %d\n", nrow(D_K), ncol(D_K)))
    cat(sprintf("  QC D_K : sym_err = %.3e | diag_err = %.3e\n", sym_err, diag_err))
  }

  list(
    D_recon = D_recon,
    D_Z = D_Z,
    sigma_f = sigma_f,
    sigma_s = sigma_s,
    K_f = K_f,
    K_s = K_s,
    K_hybride = K_hybride,
    D_K = D_K,
    sym_err = sym_err,
    diag_err = diag_err
  )
}

# Exécution standard si script sourcé
res_03b <- run_hybrid_kernel_distances(verbose = TRUE)
D_recon <- res_03b$D_recon
D_Z <- res_03b$D_Z
K_hybride <- res_03b$K_hybride
D_K <- res_03b$D_K

# ----------------------------------------------------------------------------
# BLOC TEST ISOLÉ (Canadian Weather, seed=42)
# ----------------------------------------------------------------------------
if (sys.nframe() == 0) {
  cat("\n=== TEST ISOLÉ 03b : Canadian Weather (seed=42) ===\n")
  source("setup.R")
  set.seed(42)
  DATASET <- "canadian"

  source("src/00_preprocess.R")
  source("src/01_lissage.R")
  source("src/02_fpca.R")
  source("src/02b_pca_hybride_reconstruction.R")

  res_test <- run_hybrid_kernel_distances(verbose = TRUE)
  D_K_test <- res_test$D_K

  pam_fit <- pam(as.dist(D_K_test), k = 4, diss = TRUE)
  ari <- adjustedRandIndex(pam_fit$clustering, as.integer(regions))

  cat(sprintf("[TEST] dim(D_K) = %d x %d\n", nrow(D_K_test), ncol(D_K_test)))
  cat(sprintf("[TEST] nombre de NA dans D_K = %d\n", sum(is.na(D_K_test))))
  cat(sprintf("[TEST] ARI final = %.6f\n", ari))
}
