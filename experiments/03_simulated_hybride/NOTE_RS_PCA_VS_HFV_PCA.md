# RS-PCA vs HFV-PCA dans ce projet

## But de cette note
Clarifier la difference conceptuelle entre RS-PCA et HFV-PCA, puis situer la
chaine implementee (`02b -> 03b`) dans ce cadre, avec la **nomenclature** :

- **`ω`** : parametre de la distance ponderee **Dw(α, ω)** (etape 03 / strategie B) —
  melange **fonctionnel vs vectoriel** dans la matrice de distances pour PAM.
- **`r`** : ratio **V_F / V_Y** dans l’ACP hybride (etape 02b) — pondération
  **fonctionnel / vectoriel entre blocs** avant fusion dans `chi`.
  Ce n’est **pas** le meme objet que **ω** (distance Dw).

## Etat actuel du code (`02b`) — aligné sur `02b_pca_hybride_reconstruction.R`

Implémentation **HFV « double ACP préalable »** :

- **η** = scores FPCA (étape 02).
- **Z_std** = `scale(Z)` ; **ACP** (`prcomp`) sur **Z_std** → scores **γ** (`gamma`), dimension **J** (variance cumulée ≥ 95 %, plafond `min(p, n−1, 50)`).
- **V_F** = tr(Cov(**η**)), **V_Y** = tr(Cov(**γ**)).
- **r** = V_F / V_Y ; **γ_pond** = √r · **γ** ; **χ** = [**η** | **γ_pond**].
- ACP de fusion sur **χ** → **ρ** ; reconstruction L² des courbes (loadings, `phi_k`, etc.).

La nomenclature **r** / **ω** reste inchangée (**r** = pondération inter-blocs dans l’espace des scores ; **ω** = mélange dans **Dw**).

## RS-PCA (route representation)
- La partie fonctionnelle est encodee dans un espace fini (coefficients de base
  ou scores).
- On concatene cette representation avec la partie vectorielle.
- On applique ensuite une PCA matricielle standard sur l’objet concatene.
- C’est simple, operationnel et interpretable, surtout quand la partie
  vectorielle est de petite dimension.

## HFV-PCA (route espace hybride)
- On pose un espace hybride fonctionnel + vectoriel avec un produit scalaire
  pondere.
- Le parametre de pondération inter-blocs HFV est explicite (**`r`** dans
  `02b`, calcule par **V_F / V_Y** avec **V_Y** sur les scores **γ** de la PCA vectorielle).
- La theorie est formulee via un operateur de covariance hybride et sa
  decomposition spectrale.
- En pratique, on passe par une troncature (FPCA/PCA) puis une matrice bloc de
  covariance.

## Position de la chaine actuelle (02b -> 03b)
- `02b_pca_hybride_reconstruction.R` :
  - recupere `eta` (scores FPCA),
  - PCA sur **Z_std** → `gamma`, **J**, **V_Y** = tr(Cov(gamma)),
  - calcule **`r = V_F / V_Y`** et construit **`chi`**,
  - ACP de fusion et reconstruction fonctionnelle.
- `03b_distances_noyaux_hybrides.R` :
  - calcule une distance L2 sur les courbes reconstruites,
  - combine un noyau fonctionnel et un noyau vectoriel sur **Z** original standardise,
  - convertit la similarite hybride en distance `D_K`.

Conclusion : la chaine est une implementation pratique **HFV** (double ACP + fusion + reconstruction), comparable en esprit a une logique RS par concatenation puis reduction, mais avec **pondération explicite** **`r`** entre blocs fonctionnel et vectoriel (scores **γ**, pas colonnes brutes de **Z_std** seules).

## Simulateur `Cas2_deriv` (optionnel)
- **`P_Z_SIM`** : dimension `p` du bloc vectoriel (defaut **20** ; ex. 30 pour
  un espace Z plus riche), voir `src/00_preprocess_simulated.R`.
- **`PCA_Z_Q`** : si defini (entier &lt; `p`), ACP classique sur `Z` avant le
  pipeline (reduction de dimension / bruit) — **distinct** de la PCA sur **Z_std** dans **02b**.
