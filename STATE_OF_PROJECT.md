# État du projet — Cerveau externe pour l'Agent

**Dernière mise à jour** : Mars 2026  
**Usage** : L'Agent doit lire ce fichier avant chaque action complexe pour réaligner sa micro-tâche avec la macro-vision du projet.

### Nomenclature des paramètres (ne pas confondre)

| Symbole | Où | Rôle |
|---------|-----|------|
| **ω** | Étape 03, stratégie B | Paramètre **Dw(α, ω)** : mélange fonctionnel / vectoriel dans la **matrice de distances**. |
| **r** | Étape 02b | Ratio HFV **V_F / V_Y** pour pondérer les **blocs** avant fusion dans l’ACP hybride (voir [`NOTE_RS_PCA_VS_HFV_PCA.md`](experiments/03_simulated_hybride/NOTE_RS_PCA_VS_HFV_PCA.md)). |

---

## 1. Hypothèse centrale

Comparer trois stratégies de clustering pour des données mixtes (fonctionnelles + vectorielles) :

- **Stratégie A** : FPCA → scores + variables Z → k-means
- **Stratégie B** : Distance pondérée D_w(α, ω) → PAM
- **Stratégie C** : Distance à noyaux D_K → PAM

**Question centrale** : Quelle stratégie et quels paramètres (α, ω) — choisis **sans** utiliser la vérité terrain — donnent les meilleures partitions au sens de métriques d’évaluation (dont l’ARI lorsque les labels sont disponibles) ?

**Règle de conduite (code & benchmarks)** : l’**ARI ne sert jamais** au réglage des hyperparamètres ; uniquement à l’**évaluation** des résultats. Pas de grille « oracle ARI » dans le pipeline ni dans le benchmark simulé.

## 2. Découverte principale (paradoxe)

Les **dérivées sont systématiquement informatives** lorsqu’on **évalue** avec la vérité terrain (analyses rapport / contre-exemples) : des α élevés maximisant l’ARI existent sur la grille, alors que le pipeline **sélectionne** (α, ω) par **silhouette**, ce qui peut donner α = 0 sur les datasets réels.

→ **Paradoxe** : la silhouette (critère sans vérité terrain) ne reflète pas la structure réelle ; l'ARI (avec vérité terrain) montre que les dérivées comptent.

**Implication** : Le critère de sélection des paramètres (α, ω) est le **verrou principal** pour les applications non labellisées.

## 3. Points de friction actuels

1. **Abstract vs tableau** : L'abstract parle de « α > 0.6 optimal » sans préciser « ARI-optimal ». Risque de confusion (le pipeline utilise silhouette-optimal).
2. **α « ARI-optimal »** : toute mention dans le rapport doit être clairement séparée du protocole opérationnel (pas utilisée pour choisir les paramètres dans le code).
3. **Expérience 01 (stabilité)** : La stabilité sélectionne-t-elle des (α, ω) plus proches de l'ARI-optimal que la silhouette ? → En cours d'analyse.
4. **Données SHOM** : Application non labellisée prévue — pas de vérité terrain pour valider.

## 4. Impasses logiques identifiées

- **RS-PCA** : Biais du dictionnaire (base de lissage) + mur des données incomplètes (SIR). → Abandonné au profit de l'espace hybride D_w.
- **Silhouette seule** : Ne suffit pas pour calibrer α, ω en contexte supervisé.
- **Données simulées** : générateur `Cas2_deriv` + benchmark `experiments/03_simulated_hybride/` (synthèse dans `results/SYNTHESE_BENCHMARK_SIMULE.md`) ; biblio / plan dans `docs/biblio/notes/RE_Lectures_ACP_hybride/`.

## 5. Structure du projet (rappel)

- **Pipeline** : `src/main.R` → 00_preprocess → 01_lissage → 02_fpca → 03_distances → 04_clustering → 05_visualisation  
  Chaîne optionnelle hybride : `02b_pca_hybride_reconstruction.R` → `03b_distances_noyaux_hybrides.R` (exp. simulée `experiments/03_simulated_hybride/`).
- **Expériences** : `experiments/01_stabilite/` (stabilité, nselectboot), `experiments/02_vote_criteres/`, etc.
- **Rapports** : `docs/rapport_synthese.tex` (monolithique, ~1120 lignes)
- **Compilation** : `make all` (depuis la racine) ou `make -C docs` pour LaTeX seul

## 6. Instructions pour l'Agent

Avant une action complexe (réécriture, refactor R, modification de formules) :

1. **Relire** : Hypothèse centrale (§1) et paradoxe (§2)
2. **Vérifier** : La modification ne contredit pas les découvertes documentées
3. **Consulter** : `docs/revue_rapport_synthese.md` pour les valeurs de référence
4. **Respecter** : `.cursor/rules/r_strict.mdc` pour les paramètres R
5. **Nomenclature** : **ω** (Dw) vs **r** (HFV) — voir tableau en tête de ce fichier  
6. **ARI** : métrique d’**évaluation** uniquement — jamais pour optimiser α, ω ou d’autres hyperparamètres dans le code ou les benchmarks

## 7. Prochaines étapes (priorités)

1. Clarifier l'abstract (contraste entre ce que montre une évaluation supervisée ponctuelle et le choix par silhouette — sans utiliser l'ARI comme tuner)
2. Finaliser l'analyse de l'expérience 01 (matrices de confusion, interprétation)
3. Compléter la bibliographie (Tibshirani, Calinski-Harabasz, Davies-Bouldin)
4. **HFV** : implémenter l’ACP hybride « double préalable » (PCA sur `Z` → `gamma`, puis `r` basé sur `tr(Cov(gamma))`) — plan Cursor `hfv_acp_double_préalable_*.plan.md` (référence [`docs/PLANS_CURSOR.md`](docs/PLANS_CURSOR.md)) ; mettre à jour `02b`, NOTE, benchmark simulé.
