<div align="center">
  <h1>Classification non-supervisée de données mixtes</h1>
  <p><strong>Application à des variables hybrides (Fonctionnelles & Vectorielles)</strong></p>

  ![R](https://img.shields.io/badge/R-%23276DC3.svg?style=for-the-badge&logo=r&logoColor=white)
  ![LaTeX](https://img.shields.io/badge/latex-%23008080.svg?style=for-the-badge&logo=latex&logoColor=white)
  ![Data Science](https://img.shields.io/badge/Data%20Science-Unsupervised%20Learning-orange?style=for-the-badge)
</div>

<br>

> **Auteur** : Pierre Chambet  
> **Contexte** : Projet de fin d'études M2 TRIED – CNAM (Laboratoire CEDRIC)

---

Ce dépôt présente les travaux de recherche et l'implémentation algorithmique d'un pipeline de **clustering (regroupement) pour données mixtes**.
La particularité de ces données est qu'un même individu est caractérisé simultanément par :
1. **Une composante fonctionnelle** : Une courbe continue projetée sur un espace de Hilbert $L^2$.
2. **Une composante vectorielle** : Un vecteur discret classique mesuré dans $\mathbb{R}^p$.

Le défi mathématique fondamental consiste à formuler une géométrie commune pour ces deux espaces sans occulter leur corrélation sous-jacente. Ce projet propose, compare et critique différentes architectures, notamment en proposant une **hybridation précoce (HFV)**.

## 🧠 Les 4 Stratégies de Clustering Comparées

Le pipeline R implémente et compare rigoureusement différentes méthodes de fusion :

### 1. Architectures d'Hybridation Tardive
Ces méthodes modélisent des distances sur les variations fonctionnelles ($D_0$, $D_1$) et vectorielles ($D_s$) isolément avant de les regrouper de manière *ad-hoc*.

- **Stratégie A (RS-PCA + $k$-means) :** Extraction des scores des courbes via ACP Fonctionnelle (FPCA). Banalisation des deux modalités par simple concaténation, suivie d'un algorithme $k$-means.
- **Stratégie B (Distance Additive $D_w$) :** Construction d'une dissimilarité pondérée normalisée combinant l'amplitude, les gradients et les vecteurs :

$$
D_w(\alpha,\omega) = \sqrt{\omega \underbrace{\bigl[(1-\alpha)\tilde{D}_0^2 + \alpha\tilde{D}_1^2\bigr]}_{D_p(\alpha)^2} + (1-\omega)\tilde{D}_s^2}
$$

  Suivie d'un Partitioning Around Medoids (PAM).
- **Stratégie C (Produit de Noyaux $D_K$) :** Exigence stricte de similarité simultanée via la multiplication de noyaux gaussiens $K_f \cdot K_s$ ; la matrice de dissimilarité associée s'écrit :

$$
D_K(i,j) = \sqrt{K(i,i) + K(j,j) - 2K(i,j)}
$$



### 2. L'Innovation : L'Hybridation Précoce (HFV)
*Principal apport théorique implémenté permettant d'outrepasser les limites de l'hybridation tardive.*
Au lieu de fusionner des dissimilarités en surface, l'algorithme "Hybride Fonctionnel-Vectoriel" (HFV) crée un espace hybride fondé sur la covariance croisée en amont de tout maillage géométrique.

*   On estime une **Matrice de Covariance Jointe** unifiant les variances marginales et les variances *croisées*.

$$
V = \begin{bmatrix} V_y & V_{yx} \\ V_{xy} & V_x \end{bmatrix}
$$

*   Le ratio de trace garantit l'équilibrage des énergies: $r = \frac{\mathrm{tr}(\mathrm{Cov}(\eta))}{\mathrm{tr}(\mathrm{Cov}(\gamma))}$.
*   La diagonalisation offre des composantes informées de l'influence mutuelle des variables, reconstruites dans $L^2$ avant clustering.

---

## 📊 Matériel Expérimental (Données Réelles & Simulations)

Le pipeline a été confronté à un double environnement d'évaluation stricte :
1. **3 Jeux de données réels publics** :
   - *Canadian Weather* (4 régions climatiques) : Températures annuelles $L^2$ $\times$ position/précipitations $\mathbb{R}^p$.
   - *Berkeley Growth* (2 sexes) : Courbe de croissance $L^2$ $\times$ mensurations terminales $\mathbb{R}^p$.
   - *Tecator* (3 niveaux de gras) : Spectre d'absorbance infrarouge $L^2$ $\times$ eau/protéines $\mathbb{R}^p$.
2. **4 Scénarios Simulés (Générateur intégré)** :
   Afin de moduler indépendamment la difficulté de l'espace mixte, le dépôt intègre un environnement de simulation de Monte-Carlo. Il génère 4 scénarios (nommés $S_1$ à $S_4$), allant d'un signal mixte très franc jusqu'à un écrasement quasi-total par le bruit vectoriel et fonctionnel. Cela permet de certifier rigoureusement les conditions géométriques où l'Hybridation Précoce surpasse l'état de l'art.

---

## 🔬 Découvertes & Diagnostic : "Le Verrou du Non-Supervisé"

Le projet soulève et démontre empiriquement un **Paradoxe de la Silhouette** en contexte non supervisé mixte.

> [!WARNING]
> La maximisation classique du score de Silhouette moyen pousse les algorithmes à fragmenter arbitrairement la physique des données au profit de "masses rondes irréalistes". Ce biais mène à une chute brutale de l'**Adjusted Rand Index (ARI)** face à la vérité terrain.

**Notre Solution : Instabilité Bootstrap (Fang & Wang)**  
Au lieu de faire confiance à un critère de densité interne biaisé, le dispositif recherche l'optimum hyper-paramétrique $(\alpha, \omega, k)$ en perturbant les échantillons $B$ fois (rééchantillonnage Bootstrap) et en traquant empiriquement la stabilité des frontières de décision produites.

---

## 🗂️ Structure du Projet

```text
CNAM/
├── src/                                  # Algorithmique et Core Pipeline
│   ├── 00_preprocess_*.R                 # Scripts de parsing selon le jeu
│   ├── 01_lissage.R                      # Lissage B-splines O(n) et GCV
│   ├── 02_fpca.R                         # ACP Fonctionnelle (K par variance cumulée ≥ 95%)
│   ├── 02b_pca_hybride_reconstruction.R  # Modèle théoriqe HFV matriciel
│   ├── 03_distances.R                    # Moteur de métriques (D0, D1, Dw, DK)
│   ├── 03b_distances_noyaux_hybrides.R   # Discrétisation spatiale et Noyau HFV
│   ├── 04_clustering.R                   # PAM / k-means, Benchmarker et évaluations (ARI / Silhouette)
│   ├── 05_visualisation.R                # Rendus et figures ggplot
│   └── main.R                            # Contrôleur d'orchestration
│
├── docs/                                 # Publications, Soutenances et Rapports
│   ├── rapport_stage.tex                 # Théorie statistique intégrale
│   ├── soutenance.tex                    # Slides de la présentation de projet
│   ├── Makefile                          # Logiciel de compilation LaTeX
│   └── (exports & generated)             # CSVs bruts et Tableaux TeX auto-générés
│
├── experiments/                          # Environnement de Recherche (R&D)
│   └── 01_instabilite/                   # Procédures lourdes de bootstrap et simulations (Fang & Wang)
│
├── data/                                 # Datasets (Growth, Tecator, Canadian Weather injectés via fda)
└── figures/                              # Graphiques synthétiques produits
```

---

## 🚀 Démarrage Rapide

Ce dépôt est conçu pour être 100% reproductible sans intervention complexe.

### 1. Préparation de l'environnement (R)
Le script `setup.R` configure vos dépendances (`fda`, `cluster`, `mclust`) pour la reproductibilité :
```R
source("setup.R")
```

### 2. Exécution du Pipeline complet
L'orchestrateur `main.R` exécute la chaine intégrale depuis la préparation vectorielle jusqu'aux rapports visuels :
```R
DATASET <- "canadian"   # Datasets dispos : "canadian", "tecator", "growth"
source("src/main.R")
```
*(Le script exportera ses résultats console et génèrera les graphiques dans le dossier `figures/` correspondant au jeu choisi).*

### 3. Compilation des Rapports Mathématiques
Vous pouvez produire les PDFs (`rapport_stage.pdf`, `soutenance.pdf`) incluant les théorèmes matriciels complets :
```bash
make latex
```

---
> 💡 *Note de reproductibilité : l'ensemble des grilles d'exploration simulées et études contrefactuelles ARI-opt sont hébergées dans `experiments` pour assurer intégrité et séparation des rôles (Core Vs Recherche).*
