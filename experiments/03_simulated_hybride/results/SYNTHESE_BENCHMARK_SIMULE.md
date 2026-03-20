# Synthese benchmark simule (4 scenarios x 50 seeds)

**Rapport détaillé** (méthodologie, données simulées, résultats, interprétation) : **[`../RAPPORT_DETAILLE_EXPERIENCE_03.md`](../RAPPORT_DETAILLE_EXPERIENCE_03.md)**.

## Protocole
- Donnees : simulateur `Cas2_deriv` (`K=3` classes, `n=300`, `N=60`, **`p=20` par defaut** (`P_Z_SIM`) pour un bloc `Z` plus riche et une ACP hybride (02b) plus informative ; voir `PCA_Z_Q` dans `src/00_preprocess_simulated.R`).
- Scenarios : `S1(Delta1=1, Delta2=1)`, `S2(1,0.5)`, `S3(0.5,1)`, `S4(0.5,0.5)`.
- Repetitions : 50 seeds/scenario.
- Metriques : ARI, silhouette, matrices de confusion (seed 42 par scenario).
- **Regle de conduite** : l'ARI (et toute verite terrain) ne sert **jamais** au choix des hyperparametres (`alpha`, `omega`, etc.) ; uniquement a l'**evaluation** des partitions obtenues. Les grilles utilisent la **silhouette** pour selectionner `alpha` / `(alpha, omega)` la ou un critere interne est requis.
- **Nomenclature** : `ω` = poids Dw(α,ω) fonctionnel/vectoriel (grille 04) ; `r` = ratio HFV V_F/V_Y (02b), exporte pour `DK_reconstruit` uniquement. Colonnes CSV `omega` / `r` alignees sur cette convention.

## Methodes comparees
- Baselines : `D0`, `D1`, `Df_silopt` (α par silhouette sur la grille), `Ds`.
- Historiques : `A`, `B_silopt` (selection **silhouette** sur `(α,ω)`), `C (DK ancien)`.
- Nouvelle chaine : `DK reconstruit` (`02b -> 03b`).

## Constat principal (ARI moyen par scenario)
Relancer `benchmark_all_methods_simulated.R` pour mettre a jour les CSV ; les classements par scenario sont dans `ranking_by_scenario.csv` et les moyennes dans `metrics_*`.

En resume qualitatif (datasets faciles vs difficiles) :
- `S1` / `S2` : les methodes qui exploitent fortement la forme fonctionnelle (`D1`, `Df_silopt`, `B_silopt`) restent en tete.
- `S3` / `S4` : le signal fonctionnel etant affaibli, les ecarts se resserrent ; `A` et `Ds` peuvent chuter.

## Lecture scientifique
1. Ajouter explicitement `D1` et `Df(alpha)` avec selection par **silhouette** permet de mesurer le gain lie a la forme sans utiliser l'ARI comme tuner.
2. Le **paradoxe silhouette / ARI** se lit en comparant les metriques **apres** coup sur les memes runs : la silhouette ne predit pas toujours l'ARI — sans pour autant utiliser l'ARI pour regler quoi que ce soit.
3. `DK reconstruit` propose une fusion HFV-PCA + noyaux ; sa performance se juge sur les memes regles (selection interne sans verite terrain pour les hyperparametres).

## Protocole de sorties (fichiers + KPI + lecture console)

Voir **[`../PROTOCOLE_SORTIES.md`](../PROTOCOLE_SORTIES.md)** : description de chaque CSV, définitions des métriques, et rapports automatiques via `emettre_rapport_sorties.R`.

## Sorties produites
- `metrics_all_runs.csv` : metriques run-level.
- `metrics_summary_by_scenario_method.csv` : moyennes, ecarts-types, medianes, IQR.
- `metrics_global_average_by_method.csv` : moyenne globale par methode.
- `ranking_by_scenario.csv` : classement ARI par scenario (evaluation uniquement).
- `confusion_*.csv` : matrices de confusion representatives (seed 42).
- `qc_simulated_scenarios.csv` : diagnostics qualite des donnees simulees (`r_theorique` = V_F/V_Y sur scores/ Z).

## Position RS-PCA vs HFV-PCA
Voir `../NOTE_RS_PCA_VS_HFV_PCA.md` pour la clarification methodologique et le
positionnement de la chaine `02b -> 03b`.
