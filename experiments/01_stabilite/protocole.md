# 01 — Choix de \(k\) par instabilité bootstrap (Fang & Wang, 2012)

**Source article** : Fang & Wang (2012), *Computational Statistics and Data Analysis*.  
**Implémentation** : `fpc::nselectboot` (package R **fpc**).

## Idée

Pour chaque couple **(\(\alpha\), \(\omega\))** de la grille, on construit la distance mixte **\(D_w(\alpha,\omega)\)** (même pipeline que le projet : lissage, FPCA, distances).  
On estime, par **bootstrap**, une **instabilité** des partitions pour chaque **\(k\)** candidat ; le **\(k\)** qui minimise l’instabilité moyenne est retenu pour ce couple.

Une grille complète en **(\(\alpha\), \(\omega\))** (ici 21×21, pas 0,05) permet de cartographier **\(k_\text{opt}(\alpha,\omega)\)** et de tracer des heatmaps d’instabilité.

## Scripts (pipeline actuel)

| Script | Rôle |
|--------|------|
| **`run_all_nselectboot.R`** | 4 datasets réels, \(B=150\), grille 21×21 ; **optionnel** ensuite le volet simulé (voir ci-dessous) |
| **`run_nselectboot.R`** | Un dataset (`DATASET` défini avant `source`) |
| **`run_nselectboot_simulated.R`** | **Enrichissement** : mêmes paramètres nselectboot sur données **simulées** (Cas2_deriv, aligné exp. 3). Sorties dans `results_simulated/`. N’altère pas `experiments/03_simulated_hybride/`. |
| **`analyse_nselectboot.R`** | Heatmaps `figures/nselectboot_heatmap_*.png` (4 jeux réels) |
| **`analyse_nselectboot_simulated.R`** | Heatmaps `figures/nselectboot_heatmap_sim_*.png` (scénarios S1–S4, seed de référence) |
| **`generate_confusion_nselectboot.R`** | Matrices `confusion_nselectboot_*.csv` (règle : premier couple avec \(k_\text{opt}=k_\text{vrai}\), sinon le plus proche) — **jeux réels uniquement** |
| **`run_all_complete.R`** | Enchaîne : nselectboot (4 réels **+** volet simulé activé) → heatmaps réels → heatmaps simulés → confusion |
| **`rapport_stabilite.tex`** | Rapport PDF (titre historique ; contenu = instabilité / nselectboot). Compilation : `make` |

### Volet simulé (complément contrôlé)

- **Objectif** : sous **\(k\)** et labels connus, évaluer le comportement de nselectboot avec le **même** pipeline (lissage, FPCA, distances) que sur les vraies données — **sans remplacer** l’évaluation empirique sur les 4 jeux réels.
- **Activation** : dans `run_all_nselectboot.R`, définir avant `source` : `RUN_NSELECTBOOT_SIMULATED <- TRUE`. Par défaut le volet simulé est **désactivé** pour limiter le coût CPU ; `run_all_complete.R` l’**active** automatiquement.
- **Paramètres** : `SEEDS_SIM` (défaut `1L` ; pour alignement complet avec le benchmark exp. 3 : `1:50`), `P_Z_SIM` (défaut 20), `REFERENCE_SEED` pour les CSV « référence » heatmaps (défaut 42). Voir `results_simulated/README.md`.
- **Limites** : le générateur Cas2_deriv n’est pas le « monde réel » ; les conclusions sur la méthode restent **conditionnelles** à ce modèle. Toujours présenter les résultats réels en parallèle dans le rapport.

## Ancien protocole « stabilité » (hors pipeline)

Un prototype ad hoc (bootstrap + ARI sur sous-échantillon 80 %) **sans article de référence** a été retiré du pipeline. Les scripts correspondants sont dans **`archive/`** (voir `archive/README.md`).

## Paramètres modifiables

Dans **`run_all_nselectboot.R`** : `B_NSELECT`, `ALPHAS_NB`, `OMEGAS_NB`, `K_RANGE`.

Pour le **volet simulé** (`run_nselectboot_simulated.R`) : `SEEDS_SIM` (défaut `1L` ; alignement benchmark exp. 3 : `1:50`), `P_Z_SIM`, `REFERENCE_SEED`. **Pilot / débogage** : réduire temporairement la grille en redéfinissant `ALPHAS_NB` et `OMEGAS_NB` avant `source()` (le coût est quasi proportionnel au nombre de couples \((\alpha,\omega)\) × `nselectboot`).

## Références

- Fang, Y., Wang, J. (2012). Selection of the number of clusters via the bootstrap method.
- Voir aussi `docs/biblio/theorie_stabilite_clustering.tex` pour le contexte théorique général (certains passages parlent encore de « stabilité » au sens large).
