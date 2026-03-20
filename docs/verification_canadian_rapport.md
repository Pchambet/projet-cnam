# Vérification Canadian Weather (seed=42) — LaTeX vs R

**Date** : 18 mars 2026  
**Source R** : `make test` + `src/03_distances.R`, `src/04_clustering.R`  
**Référence** : `docs/revue_rapport_synthese.md`, `docs/rapport_synthese.tex`

---

## Métriques R (Canadian Weather, seed=42)

| Métrique | Valeur R |
|----------|----------|
| D0 (Sil / ARI) | 0.285 / 0.229 |
| Ds (Sil / ARI) | 0.448 / 0.616 |
| A (Sil / ARI) | 0.395 / 0.748 |
| B silhouette-optimal (α, ω) | (0, 0) |
| B silhouette-optimal (Sil / ARI) | 0.448 / 0.616 |
| B ARI-optimal (α, ω) | (0.6, 0.8) |
| B ARI-optimal (Sil / ARI) | 0.309 / 0.898 |
| C (α) | 0 |
| C (Sil / ARI) | 0.329 / 0.686 |
| D0 rang (hors diagonale) | [3.91, 532.1] |
| D1 rang (hors diagonale) | [0.52, 6.67] |

---

## Tableau : contradictions et valeurs manquantes

| Type | Emplacement LaTeX | Valeur LaTeX | Valeur R | Remarque |
|------|-------------------|--------------|----------|----------|
| **Contradiction** | Tableau paradoxe (l.844) | Silhouette = 0.312 pour (α=0.6, ω=0.8) | 0.309 | Écart de 0.003 ; arrondi ou ancienne version |
| **Manquant** | Tableau principal (l.756) | — | (α=0, ω=0) pour B | Le tableau ne précise pas les (α, ω) par dataset ; la légende indique seulement « sélectionné par silhouette » |
| **Manquant** | Tableau principal (l.756) | — | α=0 pour C | Idem : α de la stratégie C non indiqué par dataset |
| **Manquant** | Tableau principal | — | ARI = 0.898 pour (α=0.6, ω=0.8) | L’ARI-optimal n’apparaît que dans la section paradoxe, pas dans le tableau de résultats |
| **Manquant** | Tableau paradoxe (l.849–850) | « 9 stations mal classées » / « 1 seule » | — | Ces effectifs ne sont pas calculés ni affichés par le pipeline R |

---

*Aucune correction appliquée. En attente de validation.*
