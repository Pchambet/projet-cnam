# Sorties — nselectboot sur données simulées (exp. 1 enrichie)

Alignement sur l’expérience 03 (`benchmark_all_methods_simulated.R`) : mêmes scénarios **S1–S4**, mêmes paramètres `NC_SIM`, `DELTA_SIM`, `SIGMA2_SIM`, `TAU2_SIM`, `P_Z_SIM` (défaut 20). Les scripts de l’exp. 3 ne sont **pas** modifiés.

## Fichiers

| Fichier | Contenu |
|--------|---------|
| `nselectboot_simulated_all_runs.csv` | Long format : `scenario`, `seed`, `k_vrai`, `alpha`, `omega`, `k_opt`, `instability_min`, `stabk` |
| `nselectboot_simulated_summary_by_run.csv` | Une ligne par (scénario, seed) : `prop_grid_k_match`, `median_k_opt`, `mode_k_opt` |
| `nselectboot_simulated_S*_seed<ref>.csv` | Grille complète pour heatmaps (référence : `REFERENCE_SEED`, défaut 42) |

## Coût CPU et échelle

- Une grille **21×21** × **B = 150** × plusieurs **seeds** est coûteuse (plusieurs heures selon machine).
- Par défaut : `SEEDS_SIM <- 1L` dans `run_nselectboot_simulated.R` ; pour reproduire le design exp. 3 : `SEEDS_SIM <- 1:50` avant `source()`.
- `run_all_nselectboot.R` : le volet simulé est **désactivé** par défaut (`RUN_NSELECTBOOT_SIMULATED <- FALSE`) pour garder un run rapide sur les 4 jeux réels. `run_all_complete.R` définit `RUN_NSELECTBOOT_SIMULATED <- TRUE` pour enchaîner le simulé après les 4 jeux.

## Figures

Générées par `analyse_nselectboot_simulated.R` : `figures/nselectboot_heatmap_sim_*.png`.
