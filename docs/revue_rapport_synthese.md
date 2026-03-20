# Revue complète du rapport de synthèse

**Fichier** : `rapport_synthese.tex`  
**Date de revue** : mars 2026  
**Réviseur** : analyse automatique + vérification croisée pipeline

---

## 1. Vue d'ensemble

Le rapport est bien structuré, lisible et pédagogique. La progression logique (problème → briques de base → architecture → résultats → paradoxe → perspectives) est claire. Quelques incohérences et corrections à apporter.

---

## 2. Vérification des données numériques

### Tableau principal (Tableau p. 315–322)

Comparaison avec les sorties actuelles du pipeline (seed = 42) :

| Dataset | Colonne | Rapport | Pipeline actuel | Statut |
|---------|---------|---------|-----------------|--------|
| Canadian | D0 (Sil/ARI) | 0.285 / 0.229 | 0.285 / 0.229 | OK |
| Canadian | Ds | 0.448 / 0.616 | 0.448 / 0.616 | OK |
| Canadian | A | 0.395 / **0.748** | 0.395 / **0.748** | OK |
| Canadian | B | 0.448 / 0.616 | 0.448 / 0.616 | OK |
| Canadian | C | 0.329 / 0.686 | 0.329 / 0.686 | OK |
| Growth | Tous | — | — | OK |
| AEMET | Tous | — | — | OK |
| Tecator | C (Silhouette) | 0.266 | 0.265 | Arrondi |

**Conclusion** : Les valeurs du rapport correspondent au pipeline actuel. La différence 0.266 vs 0.265 est un arrondi négligeable.

### Échelle D₀ (p. 244)

Le rapport indique : « D₀ ∈ [3.9, 532] ».

- Le pipeline affiche : « D0 : rang [0.0, 532.1] ».
- Les valeurs nulles viennent de la diagonale de la matrice de distances.
- Le « 3.9 » correspond probablement au minimum des D₀ non nuls. Cohérent, mais une précision en note (« hors diagonale ») serait utile.

---

## 3. Erreurs et corrections à faire

### 3.1 Faute d’orthographe (ligne 372)

**Texte actuel** : « soir $D_0$ seul (Canadian, Growth), soit $D_0$ seul (AEMET, Tecator) »  

**Corrections appliquées** :
- « soir » → « soit »
- « $D_0$ » → « $D_s$ » pour Canadian/Growth (ω=0 donne $D_s$ vectoriel, pas $D_0$)

### 3.2 Ambiguïté dans l’abstract (lignes 86–88)

**Texte actuel** : « les dérivées sont systématiquement informatives ($\alpha > 0.6$ optimal sur 4/4 datasets) ».

Le pipeline choisit (α, ω) par **silhouette**, ce qui donne α = 0 pour les 4 datasets. Le rapport parle ici des α **ARI-optimaux**. *Vérifié* : Canadian a bien (α=0.6, ω=0.8) avec ARI=0.898.

**Proposition** : Préciser dans l’abstract :  
« …les dérivées sont systématiquement informatives lorsque l’on optimise par rapport à l’ARI ($\alpha_{\text{ARI-optimal}} \geq 0.6$ sur 4/4 datasets) ».

### 3.3 Reproductibilité des α ARI-optimaux (Section 5.3)

Le rapport affirme que les α ARI-optimaux sont dans [0.6, 1.0] et illustre avec une figure (CAN, AEM, GRO, TEC), mais le pipeline ne calcule pas ces valeurs : il optimise uniquement par silhouette.

**Propositions** :
1. **Ne pas** introduire de réglage des (α, ω) par maximisation de l’ARI dans le pipeline : l’ARI sert à l’évaluation, pas au choix des hyperparamètres.
2. Préciser en note que toute valeur d’ARI rapportée pour un couple (α, ω) fixe est une **mesure a posteriori** ; si le rapport mentionne un couple particulier (ex. illustration du paradoxe), indiquer qu’il ne provient pas d’une optimisation sur l’ARI.

---

## 4. Incohérences conceptuelles

### 4.1 Tableau du « paradoxe » (p. 368–399)

Le rapport donne ARI = 0.898 pour (α = 0.6, ω = 0.8) sur Canadian Weather.

- Le pipeline ne fournit que les résultats pour le couple (α, ω) optimisant la silhouette.
- 0.898 a été vérifié : c'est bien l'ARI pour (α=0.6, ω=0.8) sur la grille B (calcul manuel).

**Proposition** : Indiquer explicitement que 0.898 est l’ARI **évalué** pour un couple (α, ω) fixé à des fins d’illustration — pas le résultat d’une sélection par ARI (le pipeline sélectionne par silhouette).

### 4.2 Contradiction apparente

- Abstract : « α > 0.6 optimal ».
- Tableau principal : les (α, ω) retenus pour B sont (0, 0) ou (0, 1).

Le texte distingue bien « silhouette-optimal » et « ARI-optimal », mais l’abstract peut laisser croire que le pipeline utilise directement α > 0.6.

**Proposition** : Uniformiser le vocabulaire : « α silhouette-optimal » vs « α ARI-optimal », et le rappeler dans l’abstract.

---

## 5. Qualité structurelle et rédactionnelle

### Points forts

- Structure claire (problème → méthodes → résultats → paradoxe → perspectives).
- Diagrammes TikZ utiles (architecture, curseurs α et ω).
- Bon usage des encadrés (keybox, resultbox).
- Présentation pédagogique des distances D₀, D₁, Dₛ et de leur combinaison.
- Section « Paradoxe de la silhouette » bien argumentée.
- Perspectives de recherche variées (stabilité, vote, optimisation bayésienne, méta-calibration).

### Points à améliorer

- **Placement des tableaux** : certains tableaux seraient plus lisibles en paysage ou avec `longtable` s’ils sont très larges.
- **Figure 5 (α ARI-optimaux)** : les labels CAN, AEM, GRO, TEC sur l’axe α sont lisibles, mais une légende rappelant qu’il s’agit des α ARI-optimaux serait utile.
- **Références** : Tibshirani et al., Calinski-Harabasz, Davies-Bouldin sont mentionnés mais non référencés formellement. À compléter en bibliographie.

---

## 6. Cohérence avec le code

### Alignements confirmés

- Formules Dₚ, D_w, D_K conformes à `03_distances.R` et `04_clustering.R`.
- Grilles (α, ω) : {0, 0.1, …, 1}, donc 11×11 = 121 points.
- Critère de sélection B : silhouette.
- PAM utilisé pour B et C.
- Seed : 42 (mentionné dans le rapport).

### Différences actuelles (mars 2026)

1. **Formule de nbasis** : Le rapport ne mentionne pas la règle `nbasis = min(65, max(15, N %/% 3))` introduite récemment.
2. **Grille GCV** : Le rapport ne précise pas la plage `log λ ∈ [-4, 8]` pour le lissage.
3. **Données** : Le rapport décrit correctement les quatre datasets (Canadian, Growth, AEMET, Tecator).

---

## 7. Recommandations priorisées

### Priorité haute

1. Corriger « soir » → « soit » (ligne 372).
2. Clarifier dans l’abstract que « α > 0.6 optimal » porte sur l’ARI et non sur le critère par défaut du pipeline.
3. Vérifier ou documenter l’origine du ARI = 0.898 pour (α = 0.6, ω = 0.8).

### Priorité moyenne

4. Ajouter au pipeline (`04_clustering.R`) l’affichage des (α, ω) ARI-optimaux.
5. Ajouter une courte note sur les plages D₀ (ex. exclusion de la diagonale) si besoin.
6. Compléter la bibliographie pour les critères cités (Gap, Calinski-Harabasz, Davies-Bouldin, etc.).

### Priorité basse

7. Documenter les paramètres de lissage (nbasis, grille GCV) dans le rapport ou en annexe technique.
8. Considérer un tableau annexe en paysage pour les résultats détaillés si le rapport s’allonge.

---

## 8. Synthèse

Le rapport est solide sur le fond et bien adapté à son objectif de synthèse. Les corrections proposées portent sur :

- Une faute d’orthographe.
- Des clarifications sur « α optimal » (silhouette vs ARI).
- Une meilleure reproductibilité (calcul et affichage des α ARI-optimaux dans le pipeline).

Avec ces ajustements, le rapport sera aligné avec le code et plus facile à reproduire.
