## Bibliographie ACP hybride, RS-PCA et clustering

Ce dossier rassemble les références bibliographiques et les notes de travail
liées au projet de classification non supervisée de données mixtes
fonctionnelles + vectorielles.

- **`articles/`**
  - `Jang_2021_ACPF_hybride.pdf` : article de Jang (2021) sur l'ACP
    de données hybrides fonctionnelles + vectorielles (HFV-PCA / ACPF hybride).
    C'est la référence principale pour les méthodes RS-PCA / HFV-PCA mentionnées
    dans le plan de stage.
  - `bioinformatics_19_9_1090.pdf` : Dudoit & Fridlyand (2003), *Bagging to improve
    the accuracy of a clustering procedure*, Bioinformatics, vol. 19, n° 9, p. 1090.
    Référence pour le bagging et l'instabilité en clustering.
  - Clest (2002) : Dudoit & Fridlyand, *A prediction-based resampling method for
    estimating the number of clusters*, Genome Biology, 3(7), research0036.
    Méthode de sélection de k par prédictibilité + calibration par nullité.
  - Fang & Wang (2012) : *Selection of the number of clusters via the bootstrap
    method*, Computational Statistics and Data Analysis, 56(3), 468–477.
    Implémenté dans `fpc::nselectboot`.
  - Wang (2010) : *Consistent selection of the number of clusters via crossvalidation*,
    Biometrika, 97(4), 893–904. Base théorique pour l'instabilité.

- **`theorie_stabilite_clustering.tex`** : document théorique détaillant les maths
  derrière Clest, bagging, Fang-Wang, Wang et notre protocole actuel.


- **`books/`**
  - `Ramsay_Silverman_Functional_data_analysis.pdf` : livre de Ramsay & Silverman,
    référence de base sur les données fonctionnelles et l'ACP fonctionnelle.
    Le chapitre sur l'analyse de données mixtes est directement lié à l'ACP
    hybride demandée par les encadrants.

- **`notes/`**
  - `RE_Lectures_ACP_hybride/` : notes et scripts associés aux lectures sur
    l'ACP hybride :
    - `Plan_de_simulations.pdf` : plan de simulation de données hybrides.
    - `simulations.R` : script R associé pour générer des données simulées.
    Ces fichiers servent de point de départ pour les expériences sur données
    hybrides simulées mentionnées dans le plan de stage.

Ce sous-dossier est purement documentaire : tout le code exécutable lié aux
expériences se trouve dans `src/` (briques génériques) et `experiments/`.

