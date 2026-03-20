# Réponses — Calibration des commandes de l'Agent

**Date** : Mars 2026  
**Contexte** : Projet Mars — Classification non-supervisée de données mixtes (FDA + vectoriel)

---

## 1. Questions pour calibrer les commandes de l'Agent

### 1.1 Ingestion et analyse (Rapports & Littérature)

#### Volume et format exact

| Type | Volume | Format | Emplacement |
|------|--------|--------|-------------|
| Rapports LaTeX | ~3 900 lignes, 5.8 Mo | `.tex` (source) + `.pdf` (compilé) | `docs/`, `docs/biblio/` |
| Notes bibliographiques | ~800 lignes | `.tex`, `.md` | `docs/biblio/notes/` |
| Articles de référence | ~10+ PDFs | PDF (non versionnés) | `docs/biblio/articles/` |
| Protocoles d'expériences | ~500 lignes | Markdown (`.md`) | `experiments/*/protocole.md` |

**Pas de base de données locale.** Tout est en fichiers plats (Markdown, LaTeX, PDF). Les PDFs des articles ne sont pas nettoyés : headers, tables des matières, numéros de page sont présents. Les notes `.tex` et `.md` sont relativement propres (texte rédigé manuellement).

#### Nettoyage du bruit

- **Rapports LaTeX** : propres (rédaction manuelle, pas d'extraction automatique).
- **Articles PDF** : non traités. Aucun pipeline d'OCR ou d'extraction de texte structuré.
- **revue_rapport_synthese.md** : contient une revue manuelle avec vérification des métriques (ARI, Silhouette) contre les sorties du pipeline.

#### Ce que signifie « analyser » précisément

D'après `docs/revue_rapport_synthese.md` et la structure du projet :

1. **Extraire des métriques** : ARI, Silhouette, p-values (si tests statistiques), tailles d'effet — et vérifier leur cohérence avec le code R.
2. **Identifier des contradictions** : ex. « α silhouette-optimal » vs « α ARI-optimal » (abstract vs tableau principal).
3. **Cartographier des concepts** : architecture des distances (D₀, D₁, Dₚ, D_w, D_K), stratégies A/B/C, critères de sélection (silhouette, stabilité, nselectboot).
4. **Vérifier la reproductibilité** : alignement rapport ↔ pipeline (seed=42, grilles α/ω, formules).

L'IA ne doit **pas** inventer de métriques ou modifier les formules pour « faire coller » les résultats.

---

### 1.2 Code R et exécution (La logique quantitative)

#### Modularité

**Pipeline principal (`src/`)** : **modulaire**.

```
src/main.R          → orchestration (source les étapes)
src/00_preprocess*.R → préparation données (un fichier par dataset)
src/01_lissage.R    → B-splines + GCV
src/02_fpca.R       → ACP fonctionnelle
src/03_distances.R  → D0, D1, Dp, Ds, Dw, DK
src/04_clustering.R → stratégies A, B, C + baselines
src/05_visualisation.R → figures
```

Chaque script définit des objets en mémoire (pas de fonctions encapsulées dans des packages). Le pipeline est **séquentiel** : chaque étape dépend des sorties de la précédente.

**Expériences (`experiments/01_stabilite/`)** : pipeline centré sur **nselectboot** (instabilité, Fang & Wang) : `run_all_nselectboot.R`, `analyse_nselectboot.R`, `generate_confusion_nselectboot.R`, orchestrés par `run_all_complete.R`. Anciens scripts « stabilité » ad hoc dans `archive/`. Pas de notebooks Rmd.

#### Tests unitaires

**Aucun test formel** (pas de `testthat`, pas de `expect_*`).

La validation actuelle repose sur :
- **Vérification manuelle** : `revue_rapport_synthese.md` compare les valeurs du rapport aux sorties du pipeline.
- **Reproductibilité** : `set.seed(42)` partout.
- **Convergence** : pas de vérification automatique (ex. k-means, PAM) ; l'utilisateur vérifie visuellement les figures et tableaux.

**Recommandation** : l'IA doit **ne pas** modifier les paramètres statistiques (α, ω, nbasis, λ, k) sans justification. Pour valider une modification, elle peut proposer une comparaison avant/après (ex. ARI, Silhouette) sur un dataset de référence (ex. Canadian Weather).

---

### 1.3 LaTeX et réécriture (La production)

#### Structure du manuscrit

**Monolithique** pour les rapports principaux :

- `rapport_synthese.tex` : ~1 120 lignes, **un seul fichier** (pas de `\input{}`).
- `rapport_distances.tex` : ~636 lignes.
- `rapport_canadian_weather.tex` : ~914 lignes.
- `rapport_stabilite.tex` (experiments/01_stabilite/) : ~468 lignes.

**Exception** : `docs/biblio/notes/section_rs_hfv_pca_module.tex` est conçu pour être inclus via `\input{}` dans un rapport principal, mais n'est pas utilisé dans `rapport_synthese.tex` actuellement.

#### Critère de validation pour « réécrire pour faire sens »

D'après la revue et les objectifs du projet :

1. **Densifier l'argumentation** : clarifier les liens entre formules (Dₚ, D_w, D_K) et le code (`03_distances.R`, `04_clustering.R`).
2. **Vulgariser** : garder le niveau pédagogique (keybox, resultbox, schémas TikZ) sans perdre la précision mathématique.
3. **Connecter R ↔ littérature** : ex. « α ARI-optimal ≥ 0.6 » doit être explicitement distingué de « α silhouette-optimal = 0 » ; les valeurs du rapport doivent correspondre aux sorties du pipeline.
4. **Ne pas inventer** : pas de métriques ou de résultats non calculés par le code.

**Règle d'or** : toute réécriture qui mentionne des nombres (ARI, Silhouette, α, ω) doit être vérifiable par exécution du pipeline avec `DATASET` et `seed` documentés.

---

## 2. Manipulations préalables obligatoires

### A. Bac à sable d'exécution

**État actuel** :
- `docs/Makefile` : compile uniquement le LaTeX (pas d'exécution R).
- `experiments/01_stabilite/Makefile` : compile `rapport_stabilite.tex` uniquement.
- `run_all_complete.R` : exécute les 5 étapes de l'expérience 01, mais ne compile pas le rapport.

**À mettre en place** : un point d'entrée unique (Makefile racine ou `run_all.sh`) qui :
1. Exécute `source("setup.R")` puis `source("src/main.R")` pour les 4 datasets.
2. Exécute `run_all_complete.R` pour l'expérience 01 (optionnel ou cible séparée).
3. Compile les rapports LaTeX (`docs/` et `experiments/01_stabilite/`).

→ **Fichier créé** : `Makefile` à la racine + `run_all.sh` (script bash de secours).

### B. Garde-fous (vérité terrain)

**État actuel** :
- `set.seed(42)` partout.
- `revue_rapport_synthese.md` documente les valeurs attendues (ex. Canadian : ARI=0.748 pour stratégie A, 0.898 pour (α=0.6, ω=0.8)).
- Pas de sous-échantillon figé avec résultat statistique connu.

**À mettre en place** :
1. Règle `.cursor/rules/r_strict.mdc` : interdire les modifications de variables de confusion, paramètres du modèle (α, ω, nbasis, λ, k) sans justification mathématique préalable.
2. (Optionnel) Créer `data/test_frozen/` avec un petit jeu de données et le résultat ARI/Sil attendu pour régression.

→ **Fichier créé** : `.cursor/rules/r_strict.mdc`.

### C. Cerveau externe (STATE_OF_PROJECT.md)

**État actuel** : aucun fichier centralisé. L'hypothèse et les points de friction sont dispersés dans README, protocoles, revue.

**À mettre en place** : `STATE_OF_PROJECT.md` contenant :
- Hypothèse centrale.
- Points de friction actuels (ex. paradoxe silhouette vs ARI).
- Impasses logiques identifiées.
- Instructions pour l'Agent : lire ce fichier avant chaque action complexe.

→ **Fichier créé** : `STATE_OF_PROJECT.md`.

---

## Synthèse des réponses

| Question | Réponse courte |
|----------|----------------|
| Volume rapports | ~4k lignes LaTeX, 5.8 Mo docs, format .tex/.md |
| Nettoyage bruit | Rapports propres ; PDFs non traités |
| « Analyser » | Extraire métriques, détecter contradictions, cartographier concepts, vérifier reproductibilité |
| Code R modulaire ? | Oui (src/ modulaire, experiments/ scripts dédiés) |
| Tests unitaires ? | Non (testthat absent) ; validation manuelle + seed |
| Structure LaTeX | Monolithique (un fichier par rapport) |
| Critère réécriture | Densifier + vulgariser + connecter R↔lit sans inventer |
