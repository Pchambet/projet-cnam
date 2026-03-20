# ============================================================================
# ÉTAPE 02b : ACP HYBRIDE + RECONSTRUCTION FONCTIONNELLE (HFV, double ACP)
# ============================================================================
#
# Chaîne :
#   1) eta = scores FPCA (étape 02)
#   2) Z_std = scale(Z) ; PCA(Z_std) -> gamma (J composantes)
#   3) V_F = trace(Cov(eta)), V_Y = trace(Cov(gamma))
#   4) r = V_F / V_Y ; gamma_pond = sqrt(r) * gamma
#   5) chi = [eta | gamma_pond] ; prcomp(chi) -> rho
#   6) Reconstruction L^2 des courbes via loadings eta_* et phi_k
#
# Nomenclature : r = pondération inter-blocs HFV (02b) ; ω = Dw(α,ω) (étape 03).
# ============================================================================

library(fda)

cat("\n--- Étape 02b : ACP hybride + reconstruction fonctionnelle ---\n")

# ----------------------------------------------------------------------------
# Utilitaires internes
# ----------------------------------------------------------------------------

.resolve_preprocess_file <- function(dataset) {
  switch(dataset,
    "canadian" = "src/00_preprocess.R",
    "tecator"  = "src/00_preprocess_tecator.R",
    "aemet"    = "src/00_preprocess_aemet.R",
    "growth"   = "src/00_preprocess_growth.R",
    stop(sprintf("Dataset inconnu '%s'.", dataset))
  )
}

.ensure_inputs <- function() {
  if (!exists("DATASET")) DATASET <<- "canadian"

  if (!exists("Z") || !exists("Y_brut")) {
    source(.resolve_preprocess_file(DATASET))
  }
  if (!exists("X_hat")) {
    source("src/01_lissage.R")
  }
  if (!exists("scores_fpca") || !exists("phi_k")) {
    source("src/02_fpca.R")
  }
}

# Variance cumulée ≥ seuil (comme FPCA étape 02), avec plafond J
.choose_J_pca <- function(pca_obj, n_obs, p_dim, seuil = 0.95, j_max_cap = 50L) {
  sdev <- pca_obj$sdev
  if (length(sdev) == 0L) stop("PCA vectorielle : aucune composante.")
  vars <- sdev^2
  tot <- sum(vars)
  if (tot <= 0 || !is.finite(tot)) stop("PCA vectorielle : variance totale nulle ou non finie.")
  prop <- vars / tot
  cumprop <- cumsum(prop)
  j95 <- which(cumprop >= seuil)[1]
  if (is.na(j95)) j95 <- length(sdev)
  j_cap <- min(p_dim, max(1L, n_obs - 1L), j_max_cap, length(sdev))
  J <- min(j95, j_cap)
  as.integer(max(1L, J))
}

# ----------------------------------------------------------------------------
# Pipeline ACP hybride
# ----------------------------------------------------------------------------

run_hybrid_pca_reconstruction <- function(M = NULL, verbose = TRUE) {
  .ensure_inputs()

  Y_vec <- as.data.frame(Z)
  eta <- as.matrix(scores_fpca)

  if (nrow(eta) != nrow(Y_vec)) {
    stop(sprintf(
      "Désalignement individus : nrow(eta) = %d mais nrow(Y_vec) = %d.",
      nrow(eta), nrow(Y_vec)
    ))
  }

  if (!is.null(rownames(eta)) && !is.null(rownames(Y_vec))) {
    shared_ids <- intersect(rownames(eta), rownames(Y_vec))
    if (length(shared_ids) == 0) {
      stop("Aucun identifiant commun entre eta et Y_vec.")
    }
    eta <- eta[shared_ids, , drop = FALSE]
    Y_vec <- Y_vec[shared_ids, , drop = FALSE]
  }

  n_obs <- nrow(eta)
  n_eta <- ncol(eta)
  L <- n_eta

  if (verbose) {
    cat(sprintf("  [1/6] eta (FPCA) : %d x %d (L = %d)\n", n_obs, n_eta, L))
  }

  Y_std <- scale(as.matrix(Y_vec), center = TRUE, scale = TRUE)
  if (any(!is.finite(Y_std))) {
    stop("Y_std contient des valeurs non finies (NA/Inf). Vérifier les colonnes constantes.")
  }

  p_dim <- ncol(Y_std)
  pca_z <- stats::prcomp(Y_std, center = FALSE, scale. = FALSE)
  J <- .choose_J_pca(pca_z, n_obs = n_obs, p_dim = p_dim, seuil = 0.95, j_max_cap = 50L)
  gamma <- pca_z$x[, seq_len(J), drop = FALSE]
  colnames(gamma) <- paste0("gam_PC", seq_len(J))
  rownames(gamma) <- rownames(eta)

  if (verbose) {
    cat(sprintf("  [2/6] PCA(Z_std) : J = %d composantes (var. cumulée ≥ 95%% avec plafond)\n", J))
  }

  V_y <- cov(eta)
  V_F <- sum(diag(V_y))

  V_g <- cov(gamma)
  V_Y <- sum(diag(V_g))

  if (!is.finite(V_F) || !is.finite(V_Y)) {
    stop("V_F ou V_Y n'est pas fini.")
  }
  if (V_Y <= 0) {
    stop(sprintf("V_Y doit être > 0. Valeur actuelle : %.6f", V_Y))
  }

  if (verbose) {
    cat(sprintf("  [3/6] V_F = %.6f | V_Y = tr(Cov(gamma)) = %.6f\n", V_F, V_Y))
  }

  r <- as.numeric(V_F / V_Y)
  if (length(r) != 1L || !is.finite(r) || r <= 0) {
    stop("r invalide (doit être scalaire, fini et > 0).")
  }

  gamma_pond <- gamma * sqrt(r)
  colnames(gamma_pond) <- colnames(gamma)

  if (verbose) {
    cat(sprintf("  [4/6] r = %.6f ; gamma_pond = sqrt(r) * gamma\n", r))
  }

  chi <- cbind(eta, gamma_pond)
  colnames(chi) <- c(paste0("eta_", colnames(eta)), colnames(gamma_pond))

  if (nrow(chi) != n_obs) {
    stop("Erreur de dimension : nrow(chi) != nombre d'individus.")
  }

  if (verbose) {
    cat(sprintf("  [5/6] chi : %d x %d (L+J = %d+%d)\n", nrow(chi), ncol(chi), L, J))
  }

  acp_hybride <- stats::prcomp(chi, center = TRUE, scale. = FALSE)
  rho <- acp_hybride$x
  loadings <- acp_hybride$rotation

  if (is.null(M)) {
    M_eff <- ncol(rho)
  } else {
    M_eff <- min(as.integer(M), ncol(rho))
  }
  if (M_eff < 1) stop("M doit être >= 1.")

  idx_eta <- grepl("^eta_", rownames(loadings))
  C_mat <- loadings[idx_eta, 1:M_eff, drop = FALSE]

  eta_recon <- rho[, 1:M_eff, drop = FALSE] %*% t(C_mat)
  colnames(eta_recon) <- colnames(eta)
  rownames(eta_recon) <- rownames(eta)

  mean_fd <- res_pca$meanfd
  harm_eval <- eval.fd(t_jours, phi_k[1:n_eta])
  X_recon_mat <- harm_eval %*% t(eta_recon)
  mean_eval <- as.vector(eval.fd(t_jours, mean_fd))
  X_recon_mat <- sweep(X_recon_mat, 1, mean_eval, FUN = "+")

  X_recon_fd <- smooth.basis(argvals = t_jours, y = X_recon_mat, fdParobj = fdPar(X_hat$basis, 2, 0))$fd

  if (verbose) {
    cat(sprintf("  [6/6] ACP fusion : rho = %d x %d | Recon : %d x %d\n",
                nrow(rho), ncol(rho), nrow(X_recon_mat), ncol(X_recon_mat)))
  }

  list(
    eta = eta,
    gamma = gamma,
    gamma_pond = gamma_pond,
    pca_Z = pca_z,
    L = L,
    J = J,
    V_y = V_y,
    V_g = V_g,
    V_F = V_F,
    V_Y = V_Y,
    r = r,
    Y_std = Y_std,
    chi = chi,
    acp_hybride = acp_hybride,
    rho = rho,
    loadings = loadings,
    C_mat = C_mat,
    eta_recon = eta_recon,
    X_recon_mat = X_recon_mat,
    X_recon_fd = X_recon_fd,
    t_jours = t_jours
  )
}

# ----------------------------------------------------------------------------
# Exécution standard : uniquement si source() depuis un autre script (pipeline)
# Si Rscript sur ce fichier seul -> pas d'exécution ici (voir bloc test).
# ----------------------------------------------------------------------------

if (sys.nframe() > 0) {
  res_hybride <- run_hybrid_pca_reconstruction(M = NULL, verbose = TRUE)

  eta <- res_hybride$eta
  gamma <- res_hybride$gamma
  L <- res_hybride$L
  J <- res_hybride$J
  V_y <- res_hybride$V_y
  V_g <- res_hybride$V_g
  V_F <- res_hybride$V_F
  V_Y <- res_hybride$V_Y
  r <- res_hybride$r
  chi <- res_hybride$chi
  acp_hybride <- res_hybride$acp_hybride
  rho <- res_hybride$rho
  loadings <- res_hybride$loadings
  X_recon_fd <- res_hybride$X_recon_fd
  X_recon_mat <- res_hybride$X_recon_mat

  cat(sprintf("  -> Résumé 02b : r = %.6f | L = %d | J = %d | dim(rho) = %d x %d\n",
              r, L, J, nrow(rho), ncol(rho)))
}

# ----------------------------------------------------------------------------
# Bloc de test isolé : données simulées uniquement (Rscript sur ce fichier)
# ----------------------------------------------------------------------------

if (sys.nframe() == 0) {
  cat("\n=== TEST ISOLÉ 02b : données simulées (seed=42) ===\n")
  source("setup.R")
  set.seed(42)
  SEED_SIM <- 42
  source("src/00_preprocess_simulated.R")
  source("src/01_lissage.R")
  source("src/02_fpca.R")

  res_test <- run_hybrid_pca_reconstruction(M = NULL, verbose = TRUE)

  cat(sprintf("\n[TEST] r = %.6f | L = %d | J = %d\n", res_test$r, res_test$L, res_test$J))
  cat(sprintf("[TEST] dim(rho) = %d x %d\n", nrow(res_test$rho), ncol(res_test$rho)))

  dim_init <- dim(eval.fd(t_jours, X_hat))
  dim_recon <- dim(res_test$X_recon_mat)
  cat(sprintf("[TEST] courbes init / recon (temps x n) = %d x %d vs %d x %d\n",
              dim_init[1], dim_init[2], dim_recon[1], dim_recon[2]))
}
