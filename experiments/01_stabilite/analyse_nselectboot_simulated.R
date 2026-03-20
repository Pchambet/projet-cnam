# ============================================================================
# Expérience 01 — Heatmaps nselectboot (données simulées, aligné exp. 3)
# ============================================================================
#
# Prérequis : run_nselectboot_simulated.R
#   Fichiers : nselectboot_simulated_S*_seed<REFERENCE_SEED>.csv dans results_simulated/
#
# Usage : source("experiments/01_stabilite/analyse_nselectboot_simulated.R")
#
# ============================================================================

results_dir <- "experiments/01_stabilite/results_simulated"
fig_dir     <- "experiments/01_stabilite/figures"

if (!exists("REFERENCE_SEED")) REFERENCE_SEED <- 42L

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  HEATMAPS NSELECTBOOT — SIMULÉ (par scénario, seed référence)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")

parse_stabk <- function(s) {
  v <- strsplit(as.character(s), ";")[[1]]
  v <- suppressWarnings(as.numeric(v))
  if (length(v) >= 6) v[2:6] else rep(NA, 5)
}

pattern <- sprintf("^nselectboot_simulated_S[0-9]+_seed%d\\.csv$", REFERENCE_SEED)
files <- list.files(results_dir, pattern = pattern, full.names = TRUE)

if (length(files) == 0) {
  cat(sprintf("  [SKIP] Aucun fichier référence (*_seed%d.csv) dans %s\n",
              REFERENCE_SEED, results_dir))
  cat("  Lancez d’abord run_nselectboot_simulated.R (avec ce seed inclus).\n\n")
} else {
  for (f in files) {
    base <- tools::file_path_sans_ext(basename(f))
    d <- gsub("^nselectboot_simulated_", "", base)

    nb <- read.csv(f)
    if (!"omega" %in% names(nb) && "r" %in% names(nb)) nb$omega <- nb$r
    alphas <- sort(unique(nb$alpha))
    omegas <- sort(unique(nb$omega))
    k_vals <- 2:6

    grille_nb <- expand.grid(alpha = alphas, omega = omegas, k = k_vals)
    grille_nb$instabilite <- NA

    for (i in seq_len(nrow(nb))) {
      instab <- parse_stabk(nb$stabk[i])
      for (kv in seq_along(k_vals)) {
        idx <- which(grille_nb$alpha == nb$alpha[i] &
                    grille_nb$omega == nb$omega[i] &
                    grille_nb$k == k_vals[kv])
        if (length(idx) > 0 && kv <= length(instab))
          grille_nb$instabilite[idx] <- instab[kv]
      }
    }

    png(file.path(fig_dir, paste0("nselectboot_heatmap_sim_", d, ".png")),
        width = 1200, height = 800)
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))

    for (kv in k_vals) {
      gr_k <- grille_nb[grille_nb$k == kv, ]
      mat <- matrix(NA, nrow = length(omegas), ncol = length(alphas))
      for (ia in seq_along(alphas)) {
        for (jo in seq_along(omegas)) {
          idx <- which(gr_k$alpha == alphas[ia] & gr_k$omega == omegas[jo])
          if (length(idx) > 0) mat[jo, ia] <- gr_k$instabilite[idx]
        }
      }

      image(alphas, omegas, t(mat),
            col = hcl.colors(50, "RdYlBu", rev = TRUE),
            main = sprintf("sim %s — k = %d (nselectboot)", d, kv),
            xlab = expression(alpha), ylab = expression(omega),
            xlim = range(alphas), ylim = range(omegas))
      vals <- mat[!is.na(mat)]
      if (length(vals) > 0 && diff(range(vals)) > 1e-10)
        contour(alphas, omegas, t(mat), add = TRUE, labcex = 0.7)

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
    cat(sprintf("  Figure : nselectboot_heatmap_sim_%s.png\n", d))
  }
}

cat("\n")
cat("══════════════════════════════════════════════════════════════════════\n")
cat("  Heatmaps simulées dans figures/ (préfixe nselectboot_heatmap_sim_)\n")
cat("══════════════════════════════════════════════════════════════════════\n\n")
