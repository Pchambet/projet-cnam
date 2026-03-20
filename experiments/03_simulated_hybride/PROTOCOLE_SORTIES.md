# Protocole de sorties — Expérience 03 (données simulées hybrides)

Ce document décrit **ce que produit le code** une fois l’exécution terminée, **où trouver les fichiers**, et **comment lire les KPI**. Les scripts affichent aussi un **rapport texte récapitulatif** dans la console (via `emettre_rapport_sorties.R`).

**Rapport scientifique complet** (protocole, générateur, méthodes, résultats chiffrés, discussion) : **[`RAPPORT_DETAILLE_EXPERIENCE_03.md`](RAPPORT_DETAILLE_EXPERIENCE_03.md)**.

---

## 1. Paramètre de simulation `p` (dimension du bloc vectoriel `Z`)

| Paramètre R   | Rôle |
|---------------|------|
| **`P_Z_SIM`** | Nombre de variables dans `Z` pour `Cas2_deriv` (`src/00_preprocess_simulated.R`). **Valeur par défaut : 20** (au lieu de 10) pour que l’**ACP sur `Z`** dans l’étape **02b** ait un espace vectoriel plus riche et des composantes `gamma` / dimension `J` plus parlantes. |

- Pour **`p = 30`** : avant de sourcer le pipeline ou le benchmark, exécuter `P_Z_SIM <- 30L` (ou `P_Z_SIM <<- 30L` depuis un script batch).
- Option **`PCA_Z_Q`** : réduction **optionnelle** de `Z` par ACP **avant** 01–04 (distincte de la PCA du bloc dans 02b).

---

## 2. Scripts principaux et livrables

### 2.1 `benchmark_all_methods_simulated.R`

**But** : comparer les méthodes de clustering (distances + stratégies A / B / C + chaîne `02b`→`03b`) sur 4 scénarios × 50 seeds.

| Fichier sortie | Contenu | KPI / lecture |
|----------------|---------|----------------|
| `results/metrics_all_runs.csv` | Une ligne par (scénario, seed, méthode). | **ARI**, **silhouette**, `alpha`, `omega`, `r` (HFV pour `DK_reconstruit`). |
| `results/metrics_summary_by_scenario_method.csv` | Agrégation par (scénario, méthode). | Moyenne, écart-type, médiane, IQR de l’ARI et de la silhouette. |
| `results/metrics_global_average_by_method.csv` | Moyenne sur tous les scénarios et seeds. | Vue **globale** : quelle méthode domine en moyenne. |
| `results/ranking_by_scenario.csv` | Classement par scénario. | **Ordre** des méthodes selon l’ARI moyen (évaluation uniquement). |
| `results/confusion_*.csv` | Matrices de confusion (seed 42). | Effectifs **cluster × classe** pour lecture qualitative. |

**Console** : bloc « **RAPPORT DE SORTIE — Benchmark simulé** » — liste des fichiers, définitions ARI / silhouette, tableau des moyennes globales triées par ARI, top 3 par scénario.

---

### 2.2 `qc_simulated_data.R`

**But** : diagnostiquer la qualité des données simulées **par scénario** (équilibre des classes, variance fonctionnelle vs vectorielle, couplage).

| Fichier sortie | Contenu | KPI / lecture |
|----------------|---------|----------------|
| `results/qc_simulated_scenarios.csv` | Une ligne par scénario S1–S4. | `p`, `pve_pc1`, `pve_cum_k`, **`VF_trace`**, **`VY_trace`**, **`r_theorique`**, `mean_abs_corr_eta_Z`, entropie des classes. |

**Console** : bloc « **RAPPORT DE SORTIE — QC** » + tableau imprimé.

---

### 2.3 `src/02b_pca_hybride_reconstruction.R` (exécution isolée)

En **`Rscript`**, seul le **test simulé** tourne ; la console affiche **`r`**, **`L`**, **`J`**, dimensions de **`rho`**, et tailles des matrices de reconstruction — indicateurs directs de l’ACP hybride HFV.

---

## 3. Définitions rapides des métriques

| Métrique | Interprétation |
|----------|----------------|
| **ARI** | Accord entre la partition trouvée et les classes simulées ; ∈ [-0.33, 1], 1 = identique. Sert à l’**évaluation**, pas au choix des hyperparamètres. |
| **Silhouette** | Cohésion / séparation sur la **distance utilisée** par la méthode ; utile sans labels, mais peut diverger de l’ARI (paradoxe documenté). |
| **`r` (02b)** | Ratio HFV \(V_F / V_Y\) avec \(V_Y = \mathrm{tr}(\mathrm{Cov}(\gamma))\) ; pondération fonctionnel / vectoriel dans `chi`. **≠** `ω` (mélange dans `Dw`). |

---

## 4. Synthèse rédigée

Pour une **narration** des résultats agrégés, voir `results/SYNTHESE_BENCHMARK_SIMULE.md` (à mettre à jour après changement de protocole ou de `P_Z_SIM`).

---

## 5. Réutilisation du module de rapport

```r
source("experiments/03_simulated_hybride/emettre_rapport_sorties.R")
# Après chargement des objets metrics_df, agg_global, rank_df :
emettre_rapport_benchmark(metrics_df, agg_global, rank_df, p_z = 20L)
```

Les fonctions sont pures **affichage** ; les données restent dans les CSV.
