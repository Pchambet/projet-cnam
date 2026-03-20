# ============================================================================
# ÉTAPE 3 : CALCUL DES MATRICES DE DISTANCES
# ============================================================================
#
# Architecture à deux niveaux :
#
#   NIVEAU 1 — Distance fonctionnelle Dp(α) :
#     D0 : L² entre courbes       = √∫ [X_i(t) - X_j(t)]² dt   (niveau)
#     D1 : L² entre dérivées      = √∫ [X'_i(t) - X'_j(t)]² dt (forme)
#     Dp(α) = √[(1-α)·D0² + α·D1²]    α ∈ [0,1]
#
#   NIVEAU 2 — Combinaison fonctionnel + vectoriel :
#     Ds : euclidienne sur Z standardisé
#     Dw(α,ω) = √[ω·Dp(α)² + (1-ω)·Ds²]          (stratégie B)
#     DK(α)   = √[K(i,i)+K(j,j)-2K(i,j)]          (stratégie C)
#               avec Kf = exp(-Dp(α)² / 2σ²)
#
# NOTE TECHNIQUE : les boucles O(n²) sont remplacées par des opérations
# matricielles vectorisées. Gain typique : ×5–×20 selon n.
# ============================================================================

library(fda)

cat("\n--- Étape 03 : Calcul des distances ---\n")

if (!exists("scores_fpca")) source("src/02_fpca.R")

# Grille fine pour l'évaluation des intégrales (somme de Riemann)
grille_fine <- seq(rangeval[1], rangeval[2], length.out = 1000)
delta <- (rangeval[2] - rangeval[1]) / (length(grille_fine) - 1)

# ═══════════════════════════════════════════════════════════════════════════
# A. DISTANCES FONCTIONNELLES
# ═══════════════════════════════════════════════════════════════════════════

# ─── D0 : distance L² entre courbes ───
# Mesure la proximité globale des profils
cat("  Calcul de D0 (L² entre courbes)...\n")
eval_X <- eval.fd(grille_fine, X_hat) # matrice (1000 x n)

# Vectorisé : pour chaque paire (i,j), D0² = sum((X_i - X_j)²) * delta
# = (||X_i||² + ||X_j||² - 2 <X_i, X_j>) * delta  [identité du parallélogramme]
ss_X <- colSums(eval_X^2) * delta # ||X_i||²  vecteur (n,)
inner_X <- crossprod(eval_X) * delta # Gram matrix <X_i, X_j>  (n x n)
D0_sq <- outer(ss_X, ss_X, "+") - 2 * inner_X
D0_sq[D0_sq < 0] <- 0 # correction numérique (valeurs ~-1e-14)
D0 <- sqrt(D0_sq)
rownames(D0) <- colnames(D0) <- noms_stations

# ─── D1 : distance L² entre dérivées ───
# Capture les différences de dynamique (vitesse de changement)
cat("  Calcul de D1 (L² entre dérivées)...\n")
eval_dX <- eval.fd(grille_fine, X_hat, Lfdobj = 1) # dérivées premières

ss_dX <- colSums(eval_dX^2) * delta
inner_dX <- crossprod(eval_dX) * delta
D1_sq <- outer(ss_dX, ss_dX, "+") - 2 * inner_dX
D1_sq[D1_sq < 0] <- 0
D1 <- sqrt(D1_sq)
rownames(D1) <- colnames(D1) <- noms_stations

# ─── Dp(α) : distance fonctionnelle combinée ───
# Dp(α) = √[(1-α)·(D0/max(D0))² + α·(D1/max(D1))²]
# Normalisation par le max pour remettre D0 et D1 sur la même échelle.
cat("  Calcul de Dp (D0 + D1 combinées par α)...\n")

D0_norm <- D0 / max(D0)
D1_norm <- D1 / max(D1)

compute_Dp <- function(D0_norm, D1_norm, alpha) {
  sqrt((1 - alpha) * D0_norm^2 + alpha * D1_norm^2)
}

# Grille de α à explorer
alphas <- seq(0, 1, by = 0.1)

# ═══════════════════════════════════════════════════════════════════════════
# B. DISTANCE VECTORIELLE (sur Z)
# ═══════════════════════════════════════════════════════════════════════════

cat("  Calcul de Ds (euclidienne sur Z standardisé)...\n")
Z_scaled <- scale(Z)
Ds <- as.matrix(dist(Z_scaled))
rownames(Ds) <- colnames(Ds) <- noms_stations

# ═══════════════════════════════════════════════════════════════════════════
# C. DISTANCES MIXTES (fonctionnel + vectoriel)
# ═══════════════════════════════════════════════════════════════════════════

# ─── Dw(α, ω) : distance mixte pondérée (stratégie B) ───
# ω ∈ [0,1] : poids fonctionnel vs vectoriel dans la distance (≠ r en 02b = ratio HFV).
cat("  Calcul de Dw (distance mixte, 2 paramètres α et ω)...\n")

Ds_norm <- Ds / max(Ds)

compute_Dw <- function(Dp, Ds_norm, omega) {
  sqrt(omega * Dp^2 + (1 - omega) * Ds_norm^2)
}

# Grilles pour le grid search 2D
omegas <- seq(0, 1, by = 0.1)

# ─── DK(α) : distance à noyaux (stratégie C, Ferreira & de Carvalho 2014) ───
# K = Kf · Ks (noyau produit)
# Kf = exp(-Dp(α)² / (2σf²))   noyau gaussien sur Dp
# Ks = exp(-Ds² / (2σs²))       noyau gaussien sur Z
# DK(i,j) = √[K(i,i) + K(j,j) - 2K(i,j)]
cat("  Calcul de DK (distance à noyaux, paramétré par α)...\n")

sigma_s <- median(Ds[Ds > 0])
Ks <- exp(-Ds^2 / (2 * sigma_s^2))

compute_DK <- function(Dp, Ks, sigma_f = NULL) {
  if (is.null(sigma_f)) sigma_f <- median(Dp[Dp > 0])
  Kf <- exp(-Dp^2 / (2 * sigma_f^2))
  K_prod <- Kf * Ks
  # DK vectorisé : DK² = K_ii + K_jj - 2*K_ij
  diag_K <- diag(K_prod)
  DK_sq <- outer(diag_K, diag_K, "+") - 2 * K_prod
  DK_sq[DK_sq < 0] <- 0
  DK <- sqrt(DK_sq)
  rownames(DK) <- colnames(DK) <- rownames(Dp)
  return(list(DK = DK, K = K_prod, sigma_f = sigma_f))
}

# ─── Résumé ───
cat(sprintf("  D0  : rang [%.1f, %.1f]\n", min(D0[D0 > 0]), max(D0)))
cat(sprintf("  D1  : rang [%.2f, %.2f]\n", min(D1[D1 > 0]), max(D1)))
cat(sprintf("  Ds  : rang [%.2f, %.2f]\n", min(Ds[Ds > 0]), max(Ds)))
cat("  -> Distances prêtes (Dp, Dw, DK paramétrés par α et ω).\n")
cat("  -> Prêt pour le clustering (étape 04).\n")
