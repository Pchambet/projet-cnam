#!/usr/bin/env bash
# ============================================================================
# run_all.sh — Script de secours si Make n'est pas disponible
# ============================================================================
#
# Usage : ./run_all.sh [pipeline|exp01|latex|all]
# Défaut : all
#
# ============================================================================

set -e
cd "$(dirname "$0")"

MODE="${1:-all}"

case "$MODE" in
  pipeline)
    echo ">>> Pipeline R (4 datasets)..."
    Rscript -e 'source("setup.R"); for (d in c("canadian","aemet","growth","tecator")) { DATASET <<- d; source("src/main.R") }'
    ;;
  exp01)
    echo ">>> Expérience 01..."
    Rscript -e 'source("experiments/01_stabilite/run_all_complete.R")'
    ;;
  latex)
    echo ">>> Compilation LaTeX..."
    make -C docs all
    make -C experiments/01_stabilite rapport_stabilite.pdf
    ;;
  all)
    echo ">>> Pipeline R..."
    Rscript -e 'source("setup.R"); for (d in c("canadian","aemet","growth","tecator")) { DATASET <<- d; source("src/main.R") }'
    echo ">>> Expérience 01..."
    Rscript -e 'source("experiments/01_stabilite/run_all_complete.R")'
    echo ">>> LaTeX..."
    make -C docs all
    make -C experiments/01_stabilite rapport_stabilite.pdf
    ;;
  *)
    echo "Usage: $0 [pipeline|exp01|latex|all]"
    exit 1
    ;;
esac

echo ">>> Terminé."
