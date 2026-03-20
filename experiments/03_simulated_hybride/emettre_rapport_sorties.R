# ============================================================================
# RAPPORTS DE SORTIE — Expérience 03 (simulé)
# ============================================================================
# À sourcer après benchmark / QC pour afficher un récapitulatif lisible.
# Usage : source("experiments/03_simulated_hybride/emettre_rapport_sorties.R")
# ============================================================================

#' Imprime un encadré texte (UTF-8).
.print_box <- function(lines, title = NULL) {
  w <- max(nchar(lines, type = "width", allowNA = TRUE), na.rm = TRUE)
  if (!is.null(title)) {
    cat("\n")
    cat(paste0("=", title, "\n"))
  }
  for (ln in lines) cat(ln, "\n", sep = "")
}

#' Rapport après benchmark_all_methods_simulated.R
emettre_rapport_benchmark <- function(metrics_df,
                                    agg_global,
                                    rank_df,
                                    out_dir = "experiments/03_simulated_hybride/results",
                                    p_z = NA_integer_) {
  lines <- c(
    "",
    paste0("Répertoire : ", normalizePath(out_dir, mustWork = FALSE)),
    "",
    "FICHIERS PRODUITS (à ouvrir pour le détail) :",
    "  • metrics_all_runs.csv — 1 ligne = 1 méthode × 1 seed × 1 scénario (ARI, silhouette, α, ω, r si dispo).",
    "  • metrics_summary_by_scenario_method.csv — Moyenne, écart-type, médiane, IQR par (scénario, méthode).",
    "  • metrics_global_average_by_method.csv — Moyenne globale ARI / silhouette par méthode (tous scénarios confondus).",
    "  • ranking_by_scenario.csv — Classement des méthodes par ARI moyen dans chaque scénario.",
    "  • confusion_<scenario>_<methode>.csv — Matrices de confusion (référence seed 42).",
    "",
    "KPI — LECTURE :",
    "  • ARI (Adjusted Rand Index) : accord partition vs classes simulées ; ∈ [-0.33, 1], 1 = parfait.",
    "  • Silhouette : cohésion / séparation sur la distance utilisée par la méthode ; plus élevé = mieux séparé (sans garantie de vérité terrain).",
    "  • Colonnes alpha / omega : grille distance Dw (stratégie B) ; r : ratio HFV 02b pour DK_reconstruit uniquement.",
    "",
    "MOYENNES GLOBALES (ARI ↓ tri décroissant) :"
  )
  .print_box(lines, "=== RAPPORT DE SORTIE — Benchmark simulé (exp. 03) ===")

  ag <- agg_global[order(-agg_global$ari), ]
  rownames(ag) <- NULL
  print(ag, row.names = FALSE, digits = 3)

  cat("\n--- Top 3 méthodes par scénario (ARI moyen) ---\n")
  for (sc in sort(unique(rank_df$scenario))) {
    sub <- rank_df[rank_df$scenario == sc, ]
    sub <- sub[order(sub$rank_ari), ]
    top <- head(sub, 3)
    cat(sprintf("  %s : %s\n", sc, paste(sprintf("%s (ARI=%.3f)", top$method, top$ari_mean), collapse = " > ")))
  }

  if (!is.na(p_z) && p_z > 0L) {
    cat(sprintf("\n--- Contexte simulation : dimension vectorielle p = P_Z_SIM = %d ---\n", p_z))
  }

  cat("\n--- Interprétation rapide ---\n")
  cat("Comparer les méthodes sur l’ARI si l’objectif est la récupération des classes simulées ;\n")
  cat("la silhouette compare la structure géométrique interne à chaque distance (paradoxe possible).\n")
  cat("Les hyperparamètres (α, ω) sont choisis par silhouette dans le pipeline, pas par ARI.\n\n")
}

#' Rapport après qc_simulated_data.R
emettre_rapport_qc <- function(qc_df,
                               path_csv = "experiments/03_simulated_hybride/results/qc_simulated_scenarios.csv") {
  lines <- c(
    "",
    paste0("Fichier : ", path_csv),
    "",
    "COLONNES (diagnostic par scénario Δ) :",
    "  • n, N, p, K — taille échantillon, longueur temporelle, dim. Z, nb classes.",
    "  • class_entropy — entropie des effectifs (équilibre des classes).",
    "  • pve_pc* — variance expliquée FPCA (structure fonctionnelle).",
    "  • VF_trace / VY_trace — tr(Cov(eta)) et tr(Cov(Z_std)) ; r_theorique = VF/VY (ordre de grandeur couplage).",
    "  • mean_abs_corr_eta_Z — couplage linéaire moyen scores FPCA ↔ Z standardisé.",
    "",
    "TABLE :"
  )
  .print_box(lines, "=== RAPPORT DE SORTIE — QC données simulées ===")
  print(qc_df, row.names = FALSE, digits = 4)
  cat("\n")
}
