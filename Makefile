# ============================================================================
# Makefile — Point d'entrée unique pour le projet Mars
# ============================================================================
#
# Usage (depuis la racine du projet) :
#   make all          # Pipeline R (4 datasets) + expérience 01 + compilation LaTeX
#   make pipeline     # Uniquement pipeline R (figures, tableaux)
#   make exp01        # Uniquement expérience 01 (stabilité + nselectboot)
#   make latex        # Uniquement compilation des rapports LaTeX
#   make test         # Vérification rapide (1 dataset + LaTeX)
#
# L'Agent Cursor doit utiliser : make all (ou make test pour un cycle court)
# En cas d'erreur, lire le log du terminal et corriger.
#
# ============================================================================

.PHONY: all pipeline exp01 latex test clean help

# --- Cibles principales ---
all: pipeline exp01 latex

# Pipeline R : 4 datasets → figures/ + sorties console
pipeline:
	@echo ">>> Pipeline R (4 datasets)..."
	@Rscript -e 'source("setup.R"); for (d in c("canadian","aemet","growth","tecator")) { DATASET <<- d; source("src/main.R") }'
	@echo ">>> Pipeline terminé."

# Expérience 01 : stabilité + nselectboot + analyses + matrices de confusion
exp01:
	@echo ">>> Expérience 01 (stabilité, nselectboot)..."
	@Rscript -e 'setwd("."); source("experiments/01_stabilite/run_all_complete.R")'
	@echo ">>> Expérience 01 terminée."

# Compilation LaTeX : docs/ + experiments/01_stabilite/
latex:
	@echo ">>> Compilation LaTeX (docs/)..."
	@$(MAKE) -C docs all
	@echo ">>> Compilation LaTeX (experiments/01_stabilite/)..."
	@$(MAKE) -C experiments/01_stabilite
	@echo ">>> LaTeX terminé."

# Test rapide : 1 dataset (canadian) + compilation rapport synthèse
test:
	@echo ">>> Test rapide (canadian + LaTeX)..."
	@Rscript -e 'source("setup.R"); DATASET <<- "canadian"; source("src/main.R")'
	@$(MAKE) -C docs rapport_synthese.pdf
	@echo ">>> Test terminé."

# Nettoyage des artefacts LaTeX
clean:
	$(MAKE) -C docs clean
	$(MAKE) -C experiments/01_stabilite clean

help:
	@echo "Cibles disponibles :"
	@echo "  make all      — Pipeline R + Exp01 + LaTeX (complet)"
	@echo "  make pipeline — Pipeline R uniquement (4 datasets)"
	@echo "  make exp01    — Expérience 01 uniquement"
	@echo "  make latex    — Compilation LaTeX uniquement"
	@echo "  make test     — Test rapide (canadian + rapport_synthese)"
	@echo "  make clean    — Nettoyage .aux .log .out .toc"
