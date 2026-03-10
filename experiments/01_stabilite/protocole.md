# 01 — Sélection par stabilité du clustering

**Source** : Rapport de synthèse, § 8.1 (Perspectives de recherche)

## Idée

Le bon triplet \((k, \alpha, \omega)\) est celui qui produit des clusters **reproductibles** lorsqu'on perturbe les données. On mesure la stabilité par bootstrap.

## Protocole (rapport)

Pour chaque triplet \((k, \alpha, \omega)\) dans la grille :

1. Calculer la partition \(\mathcal{C}\) sur les \(n\) individus complets (PAM avec \(\Dw(\alpha, \omega)\)).
2. Répéter \(B\) fois (\(B \approx 50\)) :
   - Sous-échantillonner 80 % des individus aléatoirement.
   - Recalculer la distance et la partition \(\mathcal{C}_b^*\) sur ce sous-échantillon.
   - Mesurer l'ARI entre \(\mathcal{C}\) (restreinte aux individus présents) et \(\mathcal{C}_b^*\).
3. Stabilité = \(\frac{1}{B}\sum_{b=1}^{B} \text{ARI}(\mathcal{C}, \mathcal{C}_b^*)\).

Le triplet avec stabilité maximale est retenu.

## Hypothèse à tester

La stabilité sélectionne-t-elle des \((\alpha, \omega)\) plus proches de l'ARI-optimal que la silhouette ?

## Grille (pour prototype rapide)

- \(\alpha, \omega \in \{0, 0.5, 1\}\) ou grille fine \(\{0, 0.1, \ldots, 1\}\)
- \(k \in \{2, 3, 4, 5, 6\}\)
- \(B = 10\) (prototype) → 50 (final)

## Critères de succès

- Sur Canadian Weather : stabilité choisit \((\alpha, \omega)\) proche de \((0.6, 0.8)\) ?
- ARI obtenu avec le triplet stabilité-optimal \(\geq\) ARI silhouette-optimal ?

## Statut

- [x] Protocole validé
- [x] Implémentation (`run_stabilite.R`)
- [x] Résultats Canadian (prototype) : stabilité-optimal (k=2, α=0.5, ω=0.5)
- [ ] Résultats sur les 4 datasets (grille fine)
