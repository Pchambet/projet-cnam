# ============================================================================
# Expérience 01 — Heatmaps nselectboot (instabilité par α, ω, k)
# ============================================================================
#
# Les heatmaps stabilité montrent : stabilité (moyenne ARI) par (α, ω) et k.
# Les heatmaps nselectboot montrent : instabilité par (α, ω) et k.
# Pour comparer : stabilité = "haut = bon", instabilité = "bas = bon".
#
# Usage : source("experiments/01_stabilite/analyse_nselectboot.R")
# Prérequis : run_all_nselectboot.R (nselectboot_*.csv dans results/)
#
# ============================================================================

results_dir <- "experiments/01_stabilite/results"
fig_dir     <- "experiments/01_stabilite/figures"
DATASETS    <- c("canadian", "aemet", "growth", "tecator")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  HEATMAPS NSELECTBOOT (instabilité)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

# Parser stabk "NA;0.03;0.08;..." → vecteur instabilité pour k=2,3,4,5,6
parse_stabk <- function(s) {
  v <- strsplit(as.character(s), ";")[[1]]
  v <- suppressWarnings(as.numeric(v))
  # Premier = k=1 (NA), indices 2-6 = k=2,3,4,5,6
  if (length(v) >= 6) v[2:6] else rep(NA, 5)
}

for (d in DATASETS) {
  f <- file.path(results_dir, paste0("nselectboot_", d, ".csv"))
  if (!file.exists(f)) {
    cat(sprintf("  [SKIP] %s : fichier absent\n", d))
    next
  }

  nb <- read.csv(f)
  if (!"omega" %in% names(nb) && "r" %in% names(nb)) nb$omega <- nb$r
  alphas <- sort(unique(nb$alpha))
  omegas <- sort(unique(nb$omega))
  k_vals <- 2:6

  # Construire grille (alpha, omega, k) → instabilité
  grille_nb <- expand.grid(alpha = alphas, omega = omegas, k = k_vals)
  grille_nb$instabilite <- NA

  for (i in 1:nrow(nb)) {
    instab <- parse_stabk(nb$stabk[i])
    for (kv in seq_along(k_vals)) {
      idx <- which(grille_nb$alpha == nb$alpha[i] &
                  grille_nb$omega == nb$omega[i] &
                  grille_nb$k == k_vals[kv])
      # parse_stabk retourne instab[1]=k=2, instab[2]=k=3, ...
      if (length(idx) > 0 && kv <= length(instab))
        grille_nb$instabilite[idx] <- instab[kv]
    }
  }

  # Heatmaps : une par k (comme stabilité)
  # Instabilité : bas = bon → palette inversée (bleu = bas, rouge = haut)
  png(file.path(fig_dir, paste0("nselectboot_heatmap_", d, ".png")),
      width = 1200, height = 800)
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

  for (kv in k_vals) {
    gr_k <- grille_nb[grille_nb$k == kv, ]
    mat <- matrix(NA, nrow = length(omegas), ncol = length(alphas))
    for (i in seq_along(alphas)) {
      for (j in seq_along(omegas)) {
        idx <- which(gr_k$alpha == alphas[i] & gr_k$omega == omegas[j])
        if (length(idx) > 0) mat[j, i] <- gr_k$instabilite[idx]
      }
    }

    image(alphas, omegas, t(mat),
          col = hcl.colors(50, "RdYlBu", rev = TRUE),  # rouge=instable, bleu=stable
          main = sprintf("%s — k = %d (nselectboot)", d, kv),
          xlab = expression(alpha), ylab = expression(omega),
          xlim = range(alphas), ylim = range(omegas))
    vals <- mat[!is.na(mat)]
    if (length(vals) > 0 && diff(range(vals)) > 1e-10)
      contour(alphas, omegas, t(mat), add = TRUE, labcex = 0.7)

    # Marquer le meilleur (instabilité min) pour ce k
    gr_k_valid <- gr_k[!is.na(gr_k$instabilite), ]
    if (nrow(gr_k_valid) > 0) {
      idx_best <- which.min(gr_k_valid$instabilite)
      ab <- gr_k_valid[idx_best, ]
      points(ab$alpha, ab$omega, pch = 4, cex = 4, lwd = 5, col = "white")
      points(ab$alpha, ab$omega, pch = 4, cex = 3, lwd = 3, col = "black")
    }
  }

  plot.new()
  legend("center",
         legend = c("Instabilité (bleu=bon, rouge=mauvais)", "× = meilleur (α,ω) pour ce k"),
         pch = c(15, 4), col = c("steelblue", "black"), pt.cex = c(2, 1.5),
         bty = "n", cex = 1.1)

  par(mfrow = c(1, 1))
  dev.off()
  cat(sprintf("  Figure : nselectboot_heatmap_%s.png\n", d))
}

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  Heatmaps nselectboot générées dans figures/\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
