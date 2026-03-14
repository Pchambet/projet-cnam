# 01 — Sélection par stabilité du clustering

**Source** : Rapport de synthèse, § 8.1 (Perspectives de recherche)

## Idée

Le bon triplet **(k, α, ω)** est celui qui produit des clusters **reproductibles** lorsqu'on perturbe les données. On mesure la stabilité par bootstrap.

## Protocole (rapport)

Pour chaque triplet **(k, α, ω)** dans la grille :

1. Calculer la partition **C** sur les **n** individus complets (PAM avec la distance mixte **D_w(α, ω)**).
2. Répéter **B** fois (par exemple **B ≈ 50**) :
   - Retirer 20 % des individus (garder 80 % au hasard, sans remise).
   - Recalculer la distance et la partition **C\_b\*** sur ce sous-échantillon.
   - Mesurer l'ARI entre **C** (restreinte aux individus présents) et **C\_b\***.
3. Stabilité :  
   `Stab(k, α, ω) = (1/B) * sum_{b=1}^B ARI(C, C_b*)`

Le triplet avec stabilité maximale est retenu.

## Détail concret du protocole

Exemple : n = 35 stations, k = 2, B = 150.

**Étape 0 (une fois)** : PAM sur les 35 stations → partition C. Chaque station a un numéro de cluster (1 ou 2).

**PAM** : minimise la somme des distances de chaque objet à son médoïde (objet réel du cluster). Chaque objet est assigné au médoïde le plus proche. Permet d'utiliser une matrice de distances arbitraire.

**Étape 1 (répétée B fois)** : À chaque fois on repart des n complets. On ne retire pas 20 % du reste à chaque fois.
1. **Sous-échantillonnage** : on tire 80 % des n au hasard (sans remise). Avec n=35, on garde 28, on ignore 7. Répétition suivante : on repart des 35, on en tire 28 autres.
2. **Sous-matrice** : on extrait les distances entre les 28 stations retenues (sous-matrice 28×28).
3. **Nouvelle partition** : PAM sur ces 28 stations → partition C_b*.
4. **Restriction** : C a 35 valeurs. On ne garde que les 28 tirées → vecteur de 28. Les deux vecteurs qu'on compare font 28 (pas 35 vs 28). Code : `C_full[idx_sub]`.
5. **ARI** : deux vecteurs de longueur 28. La position j dans les deux = la même station. Comparer = pour chaque paire (i,j), vérifier si les deux partitions sont d'accord : i et j ensemble dans les deux, ou séparés dans les deux → accord ; sinon → désaccord. L'ARI agrège ces accords.

**Étape 2** : moyenne des 150 ARI → stabilité.

## Qu'est-ce que l'ARI entre deux partitions ?

**Même ensemble** : les deux partitions portent sur les mêmes m individus (ici les 28). Pour chaque individu, on a deux labels (un par partition). Les vecteurs sont alignés : position j = même individu.

**Comparer, concrètement** : pour chaque paire (i,j) parmi les m :
- **Accord** : i et j ensemble dans les deux, ou séparés dans les deux.
- **Désaccord** : ensemble dans l'une, séparés dans l'autre.

L'ARI compte les accords et désaccords, normalise et retourne un score [-1, 1]. ARI = 1 = accord parfait.

**Pourquoi ?** Structure robuste → retirer 20 % ne bouleverse pas les 28 restants → ARI élevé. Structure fragile → partition très différente → ARI faible.

## Pourquoi 80 % ? Réponse à l'objection

Objection : « 80 % change tellement les données que la partition n'a plus rien à voir. »

Réponse : C'est exactement ce qu'on teste. Si la structure est **vraiment** présente, retirer 20 % ne devrait pas inverser les regroupements des 80 % restants. ARI élevé = reproductible. ARI faible = structure fragile ou artificielle. 80 % est un compromis standard (90 % perturbe trop peu, 50 % trop).

## Hypothèse à tester

La stabilité sélectionne-t-elle des couples **(α, ω)** plus proches de l'ARI-optimal que la silhouette ?

## Grille (pour prototype rapide)

- **α, ω** ∈ {0, 0.5, 1} ou grille fine {0, 0.1, …, 1}
- **k** ∈ {2, 3, 4, 5, 6}
- **B** = 10 (prototype) → 50 (final)

## Critères de succès

- Sur Canadian Weather : la stabilité choisit-elle un couple (α, ω) proche de (0.6, 0.8) ?
- ARI obtenu avec le triplet stabilité-optimal ≥ ARI silhouette-optimal ?

## Structure des scripts

- **`run_stabilite.R`** : exécute la stabilité pour un dataset (DATASET défini avant `source()`).
- **`run_all_stability.R`** : lance les 4 datasets avec les paramètres globaux.
- **`run_nselectboot.R`** : sélection de k par Fang & Wang (2012), `fpc::nselectboot`. Minimise l'instabilité.
- **`run_clest_prototype.R`** : prototype Clest (Dudoit & Fridlyand 2002). Sans calibration H0 pour l'instant.
- **`analyse_stabilite.R`** : heatmaps **stabilité** (stabilite_heatmap_*.png).
- **`analyse_nselectboot.R`** : heatmaps **nselectboot** (nselectboot_heatmap_*.png). Instabilité par (α, ω) et k.
- **`generate_confusion.R`** : matrices de confusion (stabilité seule).
- **`generate_confusion_both.R`** : matrices de confusion pour **les deux** méthodes (confusion_stabilite_*.csv, confusion_nselectboot_*.csv).
- **`rapport_stabilite.tex`** : rapport LaTeX résumant protocole, résultats, matrices de confusion et interprétations. Compilation : `make` depuis ce dossier.

Paramètres modifiables dans `run_all_stability.R` : **ALPHAS, OMEGAS, K_VALUES, B, SUBSAMPLE_FRAC**.

## Théorie et références

Voir **`docs/biblio/theorie_stabilite_clustering.tex`** pour les détails mathématiques.

- **Clest** (2002) : calibration par nullité ; $d_k = t_k - E[t_k | H_0]$.
- **Fang & Wang** (2012) : minimise l'instabilité (bootstrap).
- **Wang** (2010) : validation croisée, consistance asymptotique.

## Choix du sous-échantillonnage (SUBSAMPLE_FRAC)

| Valeur | Effet | Usage typique |
|--------|-------|---------------|
| **80 %** | Compromis standard : perturbation suffisante pour tester la stabilité, assez de données pour garder la structure | Recommandé par défaut |
| 90 % | Moins de perturbation, ARI plus élevés en moyenne, moins discriminatif entre triplets | Données très bruitées ou petits échantillons |
| 63 % | ≈ 1−1/e (fraction attendue dans un bootstrap) ; plus de variabilité | Certaines méthodes de stabilité |
| 50 % | Valeur par défaut de `clusterboot` (fpc) ; forte perturbation | Détecter l'instabilité, mais peut sur-pénaliser les grands k quand n est petit |

Avec n petit (ex. Canadian n=35), garder au moins ~80 % évite des sous-échantillons trop petits pour PAM à k=4 ou 6.

## Statut

- [x] Protocole validé
- [x] Implémentation (`run_stabilite.R`)
- [x] Exécution multi-datasets (`run_all_stability.R`)
- [x] Résultats sur les 4 datasets (grille α,ω ∈ {0, 0.25, 0.5, 0.75, 1}, B=15)
- [x] Analyse et heatmaps (`analyse_stabilite.R`)
- [x] Grille fine (α, ω ∈ {0, 0.1, …, 1}, 11×11) + B=50
