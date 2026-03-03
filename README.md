# Classification non-supervisée de données mixtes — Application en océanographie

> **Auteur** : Pierre  
> **Stage** : Mars 2026  
> **Encadrement** : [à compléter]

## Objectif

Comparer trois stratégies de clustering pour des données mixtes (fonctionnelles + vectorielles) en océanographie :
- **Stratégie A** : FPCA → scores + variables Z → k-means
- **Stratégie B** : Distance pondérée D_ω (fonctionnel + vectoriel) → PAM
- **Stratégie C** : Distance à noyaux D_K → PAM

## Structure

```
src/
├── 00_preprocess.R      Chargement et filtrage des données
├── 01_lissage.R         Lissage B-splines pénalisées (GCV)
├── 02_fpca.R            ACP fonctionnelle
├── 03_distances.R       D₀, D₁, D_ω, D_K
├── 04_clustering.R      PAM, évaluation (ARI, silhouette)
├── 05_visualisation.R   Figures et cartes
└── main.R               Pipeline complet
```

## Reproduire les résultats

```r
source("src/main.R")
```

## Références

- Ramsay, J.O. & Silverman, B.W. (2005). *Functional Data Analysis*. Springer.
- Ferreira, L. & de Carvalho, F. (2014). Distance-based clustering of mixed data.
