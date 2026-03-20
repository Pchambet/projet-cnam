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
| **`run_all_nselectboot.R`** | 4 datasets, \(B=150\), grille 21×21 |
| **`run_nselectboot.R`** | Un dataset (`DATASET` défini avant `source`) |
| **`analyse_nselectboot.R`** | Heatmaps `figures/nselectboot_heatmap_*.png` |
| **`generate_confusion_nselectboot.R`** | Matrices `confusion_nselectboot_*.csv` (règle : premier couple avec \(k_\text{opt}=k_\text{vrai}\), sinon le plus proche) |
| **`run_all_complete.R`** | Enchaîne : nselectboot → heatmaps → confusion |
| **`rapport_stabilite.tex`** | Rapport PDF (titre historique ; contenu = instabilité / nselectboot). Compilation : `make` |

## Ancien protocole « stabilité » (hors pipeline)

Un prototype ad hoc (bootstrap + ARI sur sous-échantillon 80 %) **sans article de référence** a été retiré du pipeline. Les scripts correspondants sont dans **`archive/`** (voir `archive/README.md`).

## Paramètres modifiables

Dans **`run_all_nselectboot.R`** : `B_NSELECT`, `ALPHAS_NB`, `OMEGAS_NB`, `K_RANGE`.

## Références

- Fang, Y., Wang, J. (2012). Selection of the number of clusters via the bootstrap method.
- Voir aussi `docs/biblio/theorie_stabilite_clustering.tex` pour le contexte théorique général (certains passages parlent encore de « stabilité » au sens large).
