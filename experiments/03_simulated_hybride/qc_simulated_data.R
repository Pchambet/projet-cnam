# ============================================================================
# QC DONNEES SIMULEES : DIAGNOSTICS DE QUALITE PAR SCENARIO
# ============================================================================

source("setup.R")

# Même défaut que 00_preprocess_simulated.R (documenté dans PROTOCOLE_SORTIES.md)
if (!exists("P_Z_SIM")) P_Z_SIM <- 20L

dir.create("experiments/03_simulated_hybride/results", recursive = TRUE, showWarnings = FALSE)

scenario_df <- data.frame(
  scenario = c("S1", "S2", "S3", "S4"),
  d1 = c(1.0, 1.0, 0.5, 0.5),
  d2 = c(1.0, 1.0, 0.5, 0.5),
  d3 = c(1.0, 0.5, 1.0, 0.5)
)

qc_rows <- list()

for (sid in seq_len(nrow(scenario_df))) {
  sc_name <- scenario_df$scenario[sid]
  DELTA_SIM <- as.numeric(scenario_df[sid, c("d1", "d2", "d3")])
  SEED_SIM <- 42
  NC_SIM <- c(100, 100, 100)
  LEN_T_SIM <- 60
  SIGMA2_SIM <- 0.2
  TAU2_SIM <- 0.2

  source("src/00_preprocess_simulated.R")
  source("src/01_lissage.R")
  source("src/02_fpca.R")

  class_tab <- table(regions)
  class_entropy <- -sum((class_tab / sum(class_tab)) * log(class_tab / sum(class_tab)))

  # Variance fonctionnelle (sans dérivée) via trace Cov(eta)
  eta <- as.matrix(scores_fpca)
  vf <- sum(diag(cov(eta)))

  # Variance vectorielle standardisée
  z_std <- scale(Z, center = TRUE, scale = TRUE)
  vy <- sum(diag(cov(z_std)))

  # Couplage hybride : moyenne des corrélations absolues scores FPCA vs Z
  corr_mat <- cor(eta, z_std)
  mean_abs_corr <- mean(abs(corr_mat))

  qc_rows[[length(qc_rows) + 1L]] <- data.frame(
    scenario = sc_name,
    delta1 = DELTA_SIM[1],
    delta2 = DELTA_SIM[3],
    n = n,
    N = N,
    p = ncol(Z),
    K = length(levels(regions)),
    class_entropy = class_entropy,
    class_min = min(class_tab),
    class_max = max(class_tab),
    pve_pc1 = res_pca$varprop[1],
    pve_pc2 = res_pca$varprop[2],
    pve_cum_k = sum(res_pca$varprop[1:K]),
    VF_trace = vf,
    VY_trace = vy,
    r_theorique = vf / vy,
    mean_abs_corr_eta_Z = mean_abs_corr
  )
}

qc_df <- do.call(rbind, qc_rows)
out_qc <- "experiments/03_simulated_hybride/results/qc_simulated_scenarios.csv"
write.csv(qc_df, file = out_qc, row.names = FALSE)

cat("QC simulé écrit dans ", out_qc, "\n", sep = "")

source("experiments/03_simulated_hybride/emettre_rapport_sorties.R")
emettre_rapport_qc(qc_df, path_csv = out_qc)
