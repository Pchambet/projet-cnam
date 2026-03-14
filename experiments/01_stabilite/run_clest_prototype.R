# ============================================================================
# Prototype Clest (Dudoit & Fridlyand 2002) pour sélection de k
# ============================================================================
#
# Principe : split train/test, cluster train, construire classifieur (affectation
# au médoïde le plus proche), prédire sur test, cluster test, comparer via
# indice FM. Calibration par données sous H0 (uniforme).
#
# Usage :
#   DATASET <- "canadian"
#   source("experiments/01_stabilite/run_clest_prototype.R")
#
# ============================================================================

library(cluster)
library(mclust)

if (!exists("DATASET")) DATASET <- "canadian"
if (!exists("B_CLEST")) B_CLEST <- 20
if (!exists("B0_CLEST")) B0_CLEST <- 20
if (!exists("TRAIN_FRAC")) TRAIN_FRAC <- 0.5
if (!exists("K_RANGE_CLEST")) K_RANGE_CLEST <- 2:6
if (!exists("PMAX")) PMAX <- 0.05
if (!exists("DMIN")) DMIN <- 0.05

# Indice Fowlkes-Mallows entre deux partitions (vecteurs de labels)
# FM = TP / sqrt((TP+FP)(TP+FN)) avec table de contingence
fowlkes_mallows <- function(c1, c2) {
  n <- length(c1)
  stopifnot(length(c2) == n)
  # Matrice de co-membership : M[i,j]=1 si i et j dans même cluster
  same1 <- outer(c1, c1, "==")
  same2 <- outer(c2, c2, "==")
  tp <- sum(same1 & same2) - n   # paires (i,i) exclues
  fp <- sum(same1 & !same2)
  fn <- sum(!same1 & same2)
  if ((tp + fp) * (tp + fn) == 0) return(0)
  tp / sqrt((tp + fp) * (tp + fn))
}

preprocess_file <- switch(DATASET,
  "canadian" = "src/00_preprocess.R",
  "tecator"  = "src/00_preprocess_tecator.R",
  "aemet"    = "src/00_preprocess_aemet.R",
  "growth"   = "src/00_preprocess_growth.R",
  "src/00_preprocess.R"
)

cat(sprintf("\n=== Clest (prototype) — %s ===\n", DATASET))
cat(sprintf("  B = %d, B0 = %d, train = %.0f%%\n", B_CLEST, B0_CLEST, TRAIN_FRAC * 100))
cat(sprintf("  k ∈ {%s}, pmax = %.2f, dmin = %.2f\n\n", paste(K_RANGE_CLEST, collapse = ","), PMAX, DMIN))

set.seed(42)

source(preprocess_file)
source("src/01_lissage.R")
source("src/02_fpca.R")
source("src/03_distances.R")

n <- ncol(Y_brut)
labels_vrai <- if (exists("regions")) as.integer(regions) else NULL
k_vrai <- if (!is.null(labels_vrai)) length(unique(labels_vrai)) else NA

# Pour ce prototype, on fixe (α, ω) au meilleur connu (ex. Canadian : 0.6, 0.8)
# ou on teste une grille réduite
alpha_clest <- 0.6
omega_clest <- 0.8

Dp <- sqrt((1 - alpha_clest) * D0_norm^2 + alpha_clest * D1_norm^2)
Dw <- sqrt(omega_clest * Dp^2 + (1 - omega_clest) * Ds_norm^2)

# Clest pour un k donné
clest_for_k <- function(Dw_mat, k, B, train_frac) {
  n <- nrow(Dw_mat)
  n_train <- max(k + 1, floor(n * train_frac))
  n_test <- n - n_train
  if (n_test < k) return(list(t_k = NA, s_k = rep(NA, B)))
  s_k <- numeric(B)
  for (b in 1:B) {
    idx_train <- sort(sample(n, size = n_train, replace = FALSE))
    idx_test <- setdiff(1:n, idx_train)
    D_train <- Dw_mat[idx_train, idx_train]
    D_test_train <- Dw_mat[idx_test, idx_train]  # distances test -> train
    pam_train <- pam(as.dist(D_train), k = k, diss = TRUE)
    # Classifier : affecter chaque point test au médoïde le plus proche
    medoid_indices <- pam_train$medoids
    medoid_global <- idx_train[medoid_indices]
    dist_to_medoids <- D_test_train[, medoid_indices]
    pred_test <- apply(dist_to_medoids, 1, which.min)
    # Cluster test
    D_test <- Dw_mat[idx_test, idx_test]
    pam_test <- pam(as.dist(D_test), k = k, diss = TRUE)
    s_k[b] <- fowlkes_mallows(pred_test, pam_test$clustering)
  }
  list(t_k = median(s_k), s_k = s_k)
}

# Clest simplifié : pas de génération H0 pour ce prototype (coûteux en calcul)
# On calcule seulement t_k et d_k = t_k - mean(t_k sous H0)
# Pour H0 rapide : permuter les lignes de la matrice de distance ? Non, ça
# ne donne pas des données uniformes. On skip H0 pour le prototype.

cat("  Calcul des t_k (sans calibration H0 pour prototype)...\n")
t_obs <- numeric(length(K_RANGE_CLEST))
names(t_obs) <- paste0("k", K_RANGE_CLEST)

for (i in seq_along(K_RANGE_CLEST)) {
  k <- K_RANGE_CLEST[i]
  res <- clest_for_k(Dw, k, B_CLEST, TRAIN_FRAC)
  t_obs[i] <- res$t_k
  cat(sprintf("    k = %d : t_k = %.4f\n", k, res$t_k))
}

# Sans H0, on prend d_k = t_k (pas de calibration)
d_obs <- t_obs
k_opt_clest <- K_RANGE_CLEST[which.max(d_obs)]
cat(sprintf("\n  k_opt (Clest, sans H0) = %d\n", k_opt_clest))
if (!is.na(k_vrai)) cat(sprintf("  k_vrai = %d\n", k_vrai))

results_dir <- "experiments/01_stabilite/results"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(data.frame(k = K_RANGE_CLEST, t_k = t_obs, d_k = d_obs),
  file.path(results_dir, sprintf("clest_prototype_%s.csv", DATASET)),
  row.names = FALSE)
