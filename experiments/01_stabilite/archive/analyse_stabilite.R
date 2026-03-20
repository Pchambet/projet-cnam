# ============================================================================
# Expérience 01 — Analyse des résultats de stabilité
# ============================================================================
#
# Usage (depuis la racine du projet) :
#   source("experiments/01_stabilite/analyse_stabilite.R")
#
# Prérequis : avoir exécuté run_all_stability.R (4 fichiers CSV dans results/)
#
# ============================================================================

results_dir <- "experiments/01_stabilite/results"
fig_dir     <- "experiments/01_stabilite/figures"
DATASETS    <- c("canadian", "aemet", "growth", "tecator")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  ANALYSE STABILITÉ — 4 datasets\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# --- 1. Charger les grilles complètes et les résumés "best" ---
grilles <- list()
bests   <- list()

for (d in DATASETS) {
  f_grille <- file.path(results_dir, paste0("stabilite_", d, ".csv"))
  f_best   <- file.path(results_dir, paste0("stabilite_", d, "_best.csv"))

  if (file.exists(f_grille)) {
    grilles[[d]] <- read.csv(f_grille)
    if (!"omega" %in% names(grilles[[d]]) && "r" %in% names(grilles[[d]])) {
      grilles[[d]]$omega <- grilles[[d]]$r
    }
    grilles[[d]]$dataset <- d
  } else {
    cat(sprintf("  [ATTENTION] Fichier manquant : %s\n", f_grille))
  }

  if (file.exists(f_best)) {
    bests[[d]] <- read.csv(f_best)
    if (!"omega" %in% names(bests[[d]]) && "r" %in% names(bests[[d]])) {
      bests[[d]]$omega <- bests[[d]]$r
    }
  } else if (d %in% names(grilles)) {
    # Fallback : extraire le meilleur de la grille (si _best absent)
    gr <- grilles[[d]]
    idx <- which.max(gr$stabilite)
    wv <- if ("omega" %in% names(gr)) gr$omega[idx] else gr$r[idx]
    bests[[d]] <- data.frame(
      dataset = d,
      k = gr$k[idx],
      alpha = gr$alpha[idx],
      omega = wv,
      stabilite = gr$stabilite[idx],
      ari = NA, n = NA, n_subsample = NA, B = NA
    )
  }
}

# Fusion des grilles et résumés
stab_all <- do.call(rbind, grilles)
best_all <- do.call(rbind, bests)
if (!"omega" %in% names(best_all) && "r" %in% names(best_all)) {
  best_all$omega <- best_all$r
}

if (is.null(stab_all)) {
  stop("Aucun fichier de résultats trouvé. Exécutez d'abord run_all_stability.R")
}

# --- 2. Tableau : meilleur triplet par dataset ---
cat("--- Meilleur triplet (stabilité max) par dataset ---\n\n")
print(best_all[, c("dataset", "k", "alpha", "omega", "stabilite", "ari")])
cat("\n")

# --- 3. Top 3 triplets par dataset ---
cat("--- Top 3 triplets par dataset ---\n\n")
for (d in DATASETS) {
  if (!d %in% names(grilles)) next
  gr <- grilles[[d]]
  gr_ord <- gr[order(-gr$stabilite), ]
  top3 <- head(gr_ord[, c("k", "alpha", "omega", "stabilite")], 3)
  cat(sprintf("  %s :\n", d))
  print(top3, row.names = FALSE)
  cat("\n")
}

# --- 4. Heatmaps stabilité (α, ω) par dataset et par k ---
# Une figure par dataset : 5 sous-plots (k=2..6)
for (d in DATASETS) {
  if (!d %in% names(grilles)) next

  gr <- grilles[[d]]
  alphas <- sort(unique(gr$alpha))
  omegas <- sort(unique(gr$omega))
  k_vals <- sort(unique(gr$k))

  png(file.path(fig_dir, paste0("stabilite_heatmap_", d, ".png")),
      width = 1200, height = 800)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

  for (kv in k_vals) {
    gr_k <- gr[gr$k == kv, ]
    mat <- matrix(NA, nrow = length(omegas), ncol = length(alphas))
    for (i in seq_along(alphas)) {
      for (j in seq_along(omegas)) {
        idx <- which(gr_k$alpha == alphas[i] & gr_k$omega == omegas[j])
        if (length(idx) > 0) mat[j, i] <- gr_k$stabilite[idx]
      }
    }

    # image() interprète z[i,j] à (x[i], y[j]) : x = lignes, y = colonnes
    # mat a lignes=omega, colonnes=alpha → on transpose pour (alpha, omega) sur (x, y)
    image(alphas, omegas, t(mat),
          col = hcl.colors(50, "YlOrRd", rev = TRUE),
          main = sprintf("%s — k = %d", d, kv),
          xlab = expression(alpha), ylab = expression(omega),
          xlim = range(alphas), ylim = range(omegas))
    contour(alphas, omegas, t(mat), add = TRUE, labcex = 0.7)

    # Marquer le meilleur pour ce k (halo blanc + croix noire pour visibilité)
    gr_k_best <- gr_k[which.max(gr_k$stabilite), ]
    if (nrow(gr_k_best) > 0) {
      points(gr_k_best$alpha, gr_k_best$omega, pch = 4, cex = 4, lwd = 5, col = "white")
      points(gr_k_best$alpha, gr_k_best$omega, pch = 4, cex = 3, lwd = 3, col = "black")
    }
  }

  plot.new()
  legend("center",
         legend = c("Stabilité (couleur)", "× = meilleur triplet"),
         pch = c(15, 4), col = c("orange", "black"), pt.cex = c(2, 1.5),
         bty = "n", cex = 1.1)

  par(mfrow = c(1, 1))
  dev.off()
  cat(sprintf("  Figure : %s\n", file.path(fig_dir, paste0("stabilite_heatmap_", d, ".png"))))
}

# --- 5. Synthèse : comparaison ARI (stabilité vs k vérité terrain) ---
cat("\n--- Synthèse ARI (triplet stabilité-optimal) ---\n\n")
n_classes <- c(canadian = 4, aemet = 4, growth = 2, tecator = 3)
for (d in DATASETS) {
  idx <- which(best_all$dataset == d)
  if (length(idx) == 0) next
  b <- best_all[idx[1], ]
  k_vrai <- n_classes[d]
  bw <- if ("omega" %in% names(b)) b$omega else b$r
  cat(sprintf("  %-12s : k_stab=%d (k_vrai=%d), α=%.2f, ω=%.2f, Stab=%.3f, ARI=%.3f\n",
              d, b$k, k_vrai, b$alpha, bw, b$stabilite, ifelse(is.na(b$ari), NA, b$ari)))
}

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  Analyse terminée. Figures dans experiments/01_stabilite/figures/\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
