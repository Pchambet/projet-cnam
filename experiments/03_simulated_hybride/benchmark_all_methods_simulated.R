# ============================================================================
# EXPÉRIENCE 03 (SIMULÉ) : BENCHMARK COMPLET DES MÉTHODES
# ============================================================================
#
# Simulation : Cas2_deriv — P_Z_SIM par défaut = 20 (voir src/00_preprocess_simulated.R).
#   Pour p = 30 : exécuter avant ce script : P_Z_SIM <<- 30L
#
# Méthodes comparées :
# - Baselines distance : D0, D1, Df(alpha)-silopt, Ds
#   (ARI uniquement en evaluation post-clustering, jamais pour choisir alpha/omega)
# - Historiques       : A, B-silopt (selection silhouette sur la grille), C (DK ancien)
# - Nouvelle chaîne   : DK reconstruit (02b + 03b) ; colonne r = ratio HFV 02b (≠ ω Dw)
#
# Design standard :
# - 4 scénarios Delta1/Delta2 : (1,1), (1,0.5), (0.5,1), (0.5,0.5)
# - 50 seeds par scénario
#
# Sorties :
# - CSV métriques run-level + agrégées
# - CSV matrices de confusion (seed 42 par scénario)
# ============================================================================

source("setup.R")

# Dimension bloc vectoriel (alignée sur 00_preprocess par défaut)
if (!exists("P_Z_SIM")) P_Z_SIM <- 20L

if (!requireNamespace("cluster", quietly = TRUE)) stop("Package 'cluster' requis.")
if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' requis.")

dir.create("experiments/03_simulated_hybride/results", recursive = TRUE, showWarnings = FALSE)

quiet_source <- function(path) {
  tf <- tempfile(fileext = ".log")
  con <- file(tf, open = "wt")
  sink(con)
  on.exit({
    sink()
    close(con)
    unlink(tf)
  }, add = TRUE)
  source(path, local = FALSE)
}

compute_pam_metrics <- function(D, labels_true, k) {
  fit <- cluster::pam(as.dist(D), k = k, diss = TRUE)
  sil <- summary(cluster::silhouette(fit$clustering, as.dist(D)))$avg.width
  ari <- mclust::adjustedRandIndex(fit$clustering, labels_true)
  list(clustering = fit$clustering, silhouette = sil, ari = ari)
}

alphas <- seq(0, 1, by = 0.1)

scenario_grid <- list(
  S1 = c(1.0, 1.0, 1.0),
  S2 = c(1.0, 1.0, 0.5),
  S3 = c(0.5, 0.5, 1.0),
  S4 = c(0.5, 0.5, 0.5)
)

SEEDS <- 1:50
SEED_REF_CONF <- 42

metrics <- list()
confusions <- list()
i_row <- 1L

for (sc_name in names(scenario_grid)) {
  DELTA_SC <- scenario_grid[[sc_name]]
  cat(sprintf("\n=== Scenario %s / delta=(%.1f, %.1f, %.1f) ===\n",
              sc_name, DELTA_SC[1], DELTA_SC[2], DELTA_SC[3]))

  for (seed in SEEDS) {
    # Paramètres simulation exposés au preprocess simulé
    SEED_SIM <- seed
    NC_SIM <- c(100, 100, 100)
    LEN_T_SIM <- 60
    DELTA_SIM <- DELTA_SC
    SIGMA2_SIM <- 0.2
    TAU2_SIM <- 0.2

    # Pipeline simulé (fonctionnel + vectoriel)
    quiet_source("src/00_preprocess_simulated.R")
    quiet_source("src/01_lissage.R")
    quiet_source("src/02_fpca.R")
    quiet_source("src/03_distances.R")
    quiet_source("src/04_clustering.R")

    labels_true <- as.integer(regions)
    k <- length(levels(regions))

    # -----------------------------------------------------------------------
    # Baselines supplémentaires demandées : D1 + Df(alpha)
    # -----------------------------------------------------------------------
    res_D1 <- compute_pam_metrics(D1, labels_true, k)

    # Df(alpha) = Dp(alpha) sur la partie fonctionnelle uniquement
    # (normalisée via D0_norm/D1_norm comme dans le pipeline historique)
    df_grid <- data.frame(alpha = alphas, silhouette = NA_real_)
    for (ia in seq_along(alphas)) {
      a <- alphas[ia]
      Df_a <- compute_Dp(D0_norm, D1_norm, alpha = a)
      m <- compute_pam_metrics(Df_a, labels_true, k)
      df_grid$silhouette[ia] <- m$silhouette
    }
    idx_df_sil <- which.max(df_grid$silhouette)

    Df_sil <- compute_Dp(D0_norm, D1_norm, alpha = df_grid$alpha[idx_df_sil])
    res_Df_sil <- compute_pam_metrics(Df_sil, labels_true, k)

    # -----------------------------------------------------------------------
    # Nouvelle chaîne 02b -> 03b (DK reconstruit)
    # -----------------------------------------------------------------------
    quiet_source("src/02b_pca_hybride_reconstruction.R")
    quiet_source("src/03b_distances_noyaux_hybrides.R")
    res_DK_recon <- compute_pam_metrics(D_K, labels_true, k)

    # -----------------------------------------------------------------------
    # Métriques consolidées de toutes les méthodes
    # -----------------------------------------------------------------------
    run_methods <- list(
      D0 = list(sil = sil_D0, ari = ari_D0, alpha = NA, omega = NA, r = NA),
      D1 = list(sil = res_D1$silhouette, ari = res_D1$ari, alpha = NA, omega = NA, r = NA),
      Df_silopt = list(sil = res_Df_sil$silhouette, ari = res_Df_sil$ari,
                       alpha = df_grid$alpha[idx_df_sil], omega = NA, r = NA),
      Ds = list(sil = sil_Ds, ari = ari_Ds, alpha = NA, omega = NA, r = NA),
      A = list(sil = resultats_A$silhouette, ari = resultats_A$ari, alpha = NA, omega = NA, r = NA),
      B_silopt = list(sil = resultats_B$silhouette, ari = resultats_B$ari,
                      alpha = resultats_B$alpha, omega = resultats_B$omega, r = NA),
      C_DK_ancien = list(sil = resultats_C$silhouette, ari = resultats_C$ari,
                         alpha = resultats_C$alpha, omega = NA, r = NA),
      DK_reconstruit = list(sil = res_DK_recon$silhouette, ari = res_DK_recon$ari,
                            alpha = NA, omega = NA, r = r)
    )

    for (mname in names(run_methods)) {
      mm <- run_methods[[mname]]
      metrics[[i_row]] <- data.frame(
        scenario = sc_name,
        delta1 = DELTA_SC[1],
        delta2 = DELTA_SC[3],
        seed = seed,
        method = mname,
        silhouette = mm$sil,
        ari = mm$ari,
        alpha = mm$alpha,
        omega = mm$omega,
        r = mm$r
      )
      i_row <- i_row + 1L
    }

    # Matrices de confusion de référence (seed 42)
    if (seed == SEED_REF_CONF) {
      confusions[[paste0(sc_name, "_A")]] <- as.data.frame.matrix(
        table(cluster = resultats_A$labels, classe = regions)
      )
      confusions[[paste0(sc_name, "_B_silopt")]] <- as.data.frame.matrix(
        table(cluster = resultats_B$labels, classe = regions)
      )
      confusions[[paste0(sc_name, "_C_DK_ancien")]] <- as.data.frame.matrix(
        table(cluster = resultats_C$labels, classe = regions)
      )
      confusions[[paste0(sc_name, "_DK_reconstruit")]] <- as.data.frame.matrix(
        table(cluster = res_DK_recon$clustering, classe = regions)
      )
      confusions[[paste0(sc_name, "_D1")]] <- as.data.frame.matrix(
        table(cluster = res_D1$clustering, classe = regions)
      )
      confusions[[paste0(sc_name, "_Df_silopt")]] <- as.data.frame.matrix(
        table(cluster = res_Df_sil$clustering, classe = regions)
      )
    }

    if (seed %% 10 == 0) {
      cat(sprintf("  - %s: seed %d/%d\n", sc_name, seed, max(SEEDS)))
    }
  }
}

metrics_df <- do.call(rbind, metrics)

agg_mean <- aggregate(cbind(ari, silhouette) ~ scenario + method, data = metrics_df, mean)
agg_sd <- aggregate(cbind(ari, silhouette) ~ scenario + method, data = metrics_df, sd)
agg_med <- aggregate(cbind(ari, silhouette) ~ scenario + method, data = metrics_df, median)
agg_iqr <- aggregate(cbind(ari, silhouette) ~ scenario + method, data = metrics_df, IQR)

agg2 <- merge(agg_mean, agg_sd, by = c("scenario", "method"), suffixes = c("_mean", "_sd"))
agg2 <- merge(agg2, agg_med, by = c("scenario", "method"))
names(agg2)[names(agg2) == "ari"] <- "ari_median"
names(agg2)[names(agg2) == "silhouette"] <- "sil_median"
agg2 <- merge(agg2, agg_iqr, by = c("scenario", "method"))
names(agg2)[names(agg2) == "ari"] <- "ari_iqr"
names(agg2)[names(agg2) == "silhouette"] <- "sil_iqr"

out_dir <- "experiments/03_simulated_hybride/results"
write.csv(metrics_df,
          file = file.path(out_dir, "metrics_all_runs.csv"),
          row.names = FALSE)
write.csv(agg2,
          file = file.path(out_dir, "metrics_summary_by_scenario_method.csv"),
          row.names = FALSE)

# Moyennes globales par méthode + classement ARI par scénario (dérivés de metrics_df)
agg_global <- aggregate(cbind(ari, silhouette) ~ method, data = metrics_df, mean)
write.csv(agg_global,
          file = file.path(out_dir, "metrics_global_average_by_method.csv"),
          row.names = FALSE)

rank_rows <- list()
for (sc in unique(metrics_df$scenario)) {
  ms <- subset(metrics_df, scenario == sc)
  a <- aggregate(ari ~ method, data = ms, mean)
  a <- a[order(-a$ari), ]
  a$silhouette_mean <- sapply(a$method, function(mm) mean(ms$silhouette[ms$method == mm]))
  a$rank_ari <- seq_len(nrow(a))
  a$scenario <- sc
  rank_rows[[length(rank_rows) + 1L]] <- a[, c("scenario", "method", "ari", "silhouette_mean", "rank_ari")]
}
rank_df <- do.call(rbind, rank_rows)
names(rank_df)[names(rank_df) == "ari"] <- "ari_mean"
write.csv(rank_df, file.path(out_dir, "ranking_by_scenario.csv"), row.names = FALSE)

for (nm in names(confusions)) {
  write.csv(confusions[[nm]],
            file = file.path(out_dir, paste0("confusion_", nm, ".csv")),
            row.names = TRUE)
}

cat("\n=== BENCHMARK TERMINE ===\n")
cat(sprintf("Runs total: %d\n", nrow(metrics_df)))
cat(sprintf("Summary rows: %d\n", nrow(agg2)))
cat(sprintf("Outputs: %s\n", out_dir))

source("experiments/03_simulated_hybride/emettre_rapport_sorties.R")
emettre_rapport_benchmark(
  metrics_df = metrics_df,
  agg_global = agg_global,
  rank_df = rank_df,
  out_dir = out_dir,
  p_z = as.integer(P_Z_SIM)
)
