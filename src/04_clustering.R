# ============================================================================
# ÉTAPE 4 : CLUSTERING — Comparaison des 3 stratégies
# ============================================================================
#
# STRATÉGIE A : FPCA → scores ξ + variables Z → k-means sur W = [ξ | Z]
#   → Tout est ramené en vecteurs, clustering classique
#
# STRATÉGIE B : Dw(α, ω) → PAM, grid search 2D
#   → α contrôle le mix D0/D1 dans Dp (niveau vs forme)
#   → ω contrôle le mix fonctionnel/vectoriel dans Dw
#   → On cherche le (α, ω) qui maximise la silhouette (l'ARI sert uniquement à l'évaluation)
#
# STRATÉGIE C : DK(α) → PAM, grid search sur α
#   → Le noyau fonctionnel utilise Dp(α) au lieu de D0
#   → σ par heuristique de la médiane, α par silhouette
#
# ÉVALUATION :
#   - Silhouette moyenne (qualité intrinsèque des clusters)
#   - ARI (Adjusted Rand Index) vs vérité terrain
#
# ============================================================================

library(cluster)   # pam, silhouette
library(mclust)    # adjustedRandIndex

cat("\n--- Étape 04 : Clustering (3 stratégies) ---\n")

if (!exists("compute_DK")) source("src/03_distances.R")

# Nombre de clusters = nombre de classes (vérité terrain)
n_clusters <- length(levels(regions))
labels_vrai <- as.integer(regions)

# ═══════════════════════════════════════════════════════════════════════════
# STRATÉGIE A : FPCA + concaténation + k-means
# ═══════════════════════════════════════════════════════════════════════════
cat("\n  ── Stratégie A : FPCA + Z → k-means ──\n")

# On concatène les scores FPCA (ξ) et les variables Z standardisées
W <- cbind(scale(scores_fpca), Z_scaled)

# k-means (avec plusieurs départs pour éviter les minima locaux)
set.seed(42)
km_A <- kmeans(W, centers = n_clusters, nstart = 25)

sil_A <- summary(silhouette(km_A$cluster, dist(W)))$avg.width
ari_A <- adjustedRandIndex(km_A$cluster, labels_vrai)
cat(sprintf("    Silhouette = %.3f | ARI = %.3f\n", sil_A, ari_A))

resultats_A <- list(
  labels = km_A$cluster,
  silhouette = sil_A,
  ari = ari_A,
  methode = "k-means sur [scores_FPCA | Z]"
)

# ═══════════════════════════════════════════════════════════════════════════
# STRATÉGIE B : Dw(α, ω) → PAM — Grid search 2D
# ═══════════════════════════════════════════════════════════════════════════
cat("\n  ── Stratégie B : Dw(α,ω) → PAM (grid search 2D) ──\n")

# Construction de la grille 2D (α, ω)
grille_B <- expand.grid(alpha = alphas, omega = omegas)
grille_B$silhouette <- NA
grille_B$ari <- NA

for (idx in 1:nrow(grille_B)) {
  a <- grille_B$alpha[idx]
  w <- grille_B$omega[idx]
  Dp_temp <- compute_Dp(D0_norm, D1_norm, alpha = a)
  Dw_temp <- compute_Dw(Dp_temp, Ds_norm, omega = w)
  pam_temp <- pam(as.dist(Dw_temp), k = n_clusters, diss = TRUE)
  grille_B$silhouette[idx] <- pam_temp$silinfo$avg.width
  grille_B$ari[idx] <- adjustedRandIndex(pam_temp$clustering, labels_vrai)
}

# Meilleur couple (α, ω) — par silhouette uniquement (jamais par ARI pour le choix)
idx_best_B <- which.max(grille_B$silhouette)
alpha_best_B <- grille_B$alpha[idx_best_B]
omega_best_B <- grille_B$omega[idx_best_B]

Dp_best <- compute_Dp(D0_norm, D1_norm, alpha = alpha_best_B)
Dw_best <- compute_Dw(Dp_best, Ds_norm, omega = omega_best_B)
pam_B <- pam(as.dist(Dw_best), k = n_clusters, diss = TRUE)

cat(sprintf("    Meilleur α = %.1f, ω = %.1f (silhouette)\n", alpha_best_B, omega_best_B))
cat(sprintf("    Silhouette = %.3f | ARI = %.3f\n",
            grille_B$silhouette[idx_best_B],
            grille_B$ari[idx_best_B]))

resultats_B <- list(
  labels = pam_B$clustering,
  silhouette = grille_B$silhouette[idx_best_B],
  ari = grille_B$ari[idx_best_B],
  alpha = alpha_best_B,
  omega = omega_best_B,
  grid = grille_B,
  methode = sprintf("PAM sur Dw (α=%.1f, ω=%.1f)", alpha_best_B, omega_best_B)
)

# ═══════════════════════════════════════════════════════════════════════════
# STRATÉGIE C : DK(α) → PAM — Grid search sur α
# ═══════════════════════════════════════════════════════════════════════════
cat("\n  ── Stratégie C : DK(α) → PAM (grid search sur α) ──\n")

grille_C <- data.frame(alpha = alphas, silhouette = NA, ari = NA)

for (idx in seq_along(alphas)) {
  a <- alphas[idx]
  Dp_temp <- compute_Dp(D0_norm, D1_norm, alpha = a)
  res_DK <- compute_DK(Dp_temp, Ks)
  pam_temp <- pam(as.dist(res_DK$DK), k = n_clusters, diss = TRUE)
  grille_C$silhouette[idx] <- pam_temp$silinfo$avg.width
  grille_C$ari[idx] <- adjustedRandIndex(pam_temp$clustering, labels_vrai)
}

# Meilleur α
idx_best_C <- which.max(grille_C$silhouette)
alpha_best_C <- grille_C$alpha[idx_best_C]
Dp_best_C <- compute_Dp(D0_norm, D1_norm, alpha = alpha_best_C)
res_DK_best <- compute_DK(Dp_best_C, Ks)
DK <- res_DK_best$DK
pam_C <- pam(as.dist(DK), k = n_clusters, diss = TRUE)
sil_C <- pam_C$silinfo$avg.width
ari_C <- adjustedRandIndex(pam_C$clustering, labels_vrai)

cat(sprintf("    Meilleur α = %.1f\n", alpha_best_C))
cat(sprintf("    Silhouette = %.3f | ARI = %.3f\n", sil_C, ari_C))

resultats_C <- list(
  labels = pam_C$clustering,
  silhouette = sil_C,
  ari = ari_C,
  alpha = alpha_best_C,
  grid = grille_C,
  methode = sprintf("PAM sur DK (α=%.1f)", alpha_best_C)
)

# ═══════════════════════════════════════════════════════════════════════════
# BASELINES pour comparaison
# ═══════════════════════════════════════════════════════════════════════════
cat("\n  ── Baselines ──\n")

# D0 seul (purement fonctionnel, pas de Z)
pam_D0 <- pam(as.dist(D0), k = n_clusters, diss = TRUE)
sil_D0 <- pam_D0$silinfo$avg.width
ari_D0 <- adjustedRandIndex(pam_D0$clustering, labels_vrai)
cat(sprintf("    D0 seul  : Silhouette = %.3f | ARI = %.3f\n", sil_D0, ari_D0))

# Ds seul (Z seul, pas de courbes)
pam_Ds <- pam(as.dist(Ds), k = n_clusters, diss = TRUE)
sil_Ds <- pam_Ds$silinfo$avg.width
ari_Ds <- adjustedRandIndex(pam_Ds$clustering, labels_vrai)
cat(sprintf("    Ds seul  : Silhouette = %.3f | ARI = %.3f\n", sil_Ds, ari_Ds))

# ═══════════════════════════════════════════════════════════════════════════
# TABLEAU COMPARATIF
# ═══════════════════════════════════════════════════════════════════════════
cat("\n  ══════════════════════════════════════════════\n")
cat("  COMPARAISON FINALE (n_clusters =", n_clusters, ")\n")
cat("  ══════════════════════════════════════════════\n")

comparaison <- data.frame(
  Strategie  = c("Baseline D0", "Baseline Ds",
                  "A (FPCA+Z)",
                  sprintf("B (Dw α=%.1f ω=%.1f)", alpha_best_B, omega_best_B),
                  sprintf("C (DK α=%.1f)", alpha_best_C)),
  Silhouette = round(c(sil_D0, sil_Ds, sil_A,
                        resultats_B$silhouette, sil_C), 3),
  ARI        = round(c(ari_D0, ari_Ds, ari_A,
                        resultats_B$ari, ari_C), 3)
)
print(comparaison, row.names = FALSE)
cat("  -> Prêt pour la visualisation (étape 05).\n")