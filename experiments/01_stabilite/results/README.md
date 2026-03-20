# Résultats — expérience 01 (instabilité / nselectboot)

## Fichiers maintenus par le pipeline actuel

- `nselectboot_*.csv` — grille 21×21, \(B=150\), un fichier par dataset.
- `confusion_nselectboot_*.csv` — matrices de confusion (règle : voir `generate_confusion_nselectboot.R`).

Les heatmaps correspondantes sont dans `../figures/` (`nselectboot_heatmap_*.png`, non versionnées si ignorées par `.gitignore`).

## Fichiers `stabilite_*` et `confusion_stabilite_*`

Ces sorties correspondaient à l’ancien protocole « stabilité » (hors pipeline depuis le recalibrage). **Ils ne sont plus régénérés.**  
Tu peux les supprimer du dossier ou les garder comme archive locale ; ils ne sont plus documentés dans le rapport PDF de l’exp. 1.

## Autres

- `clest_prototype_canadian.csv` : prototype séparé (non inclus dans `run_all_complete.R`).
