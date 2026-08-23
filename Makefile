# Build everything for the ELP gelation paper: figures from data, then PDFs.
#
#   make env        create the conda environment (once)
#   make figures    regenerate every data-derived figure in analysis/
#   make papers     compile main.pdf, SI.pdf, cover_letter.pdf into manuscript/examples/
#   make all        figures + papers
#   make check      report undefined citations / references / missing figures
#   make clean      remove LaTeX aux files (keeps the PDFs)
#
# Figure 1 (Setup.pdf) is the approved illustrator schematic and is NOT
# generated here. Every other figure is produced from the experimental data by
# the scripts under analysis/, and is named exactly as main.tex / si.tex
# reference it -- only the directory differs, which TEXINPUTS below resolves.

SHELL := /bin/bash
ROOT  := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))

MS    := $(ROOT)/manuscript
ANA   := $(ROOT)/analysis
BUILD := $(MS)/build
OUT   := $(MS)/examples

CONDA_ENV := env-2021-LLDP
PY        := conda run --no-capture-output -n $(CONDA_ENV) python

# achemso (ACS/Nano Letters class) and captdef are not in the system TeX Live,
# so they live in a repo-local texmf tree -- this keeps the build self-contained.
export TEXMFHOME := $(MS)/texmf
# Bare figure names in the .tex resolve through here; that is what lets the
# manuscript keep short names while the files live under analysis/.
export TEXINPUTS := .:$(MS)/figures:$(ANA)/SAXS_12id:$(ANA)/SAXPCS_8id:$(ANA)/Rad_Dam_Check:
# $(BUILD) must be on BIBINPUTS: achemso writes a control file (acs-main.bib)
# next to the .aux, and BibTeX aborts all databases if it cannot open it.
export BIBINPUTS := .:$(MS):$(BUILD):
export MPLBACKEND := Agg

# --- data-derived figures: script -> outputs --------------------------------
F_12ID    := $(ANA)/SAXS_12id/Figure2_SAXS_WAXS.pdf
F_GUINIER := $(ANA)/SAXS_12id/FigureS5_Guinier.pdf
F_XPCS    := $(ANA)/SAXPCS_8id/Figure3_Isothermal_SAXPCS.pdf \
             $(ANA)/SAXPCS_8id/FigureS7_Fit_Parameters.pdf \
             $(ANA)/SAXPCS_8id/FigureS3_Calibration.pdf
F_EVOL    := $(ANA)/SAXPCS_8id/FigureS6_SAXS_Evolution.pdf
F_GRID    := $(ANA)/SAXPCS_8id/FigureS8_g2_Grid.pdf
F_RAD     := $(ANA)/Rad_Dam_Check/FigureS4_Flux_Control.pdf

FIGURES := $(F_12ID) $(F_GUINIER) $(F_XPCS) $(F_EVOL) $(F_GRID) $(F_RAD)

XPCS_DATA := $(wildcard $(ANA)/SAXPCS_8id/data/*.hdf)
CSV_DATA  := $(wildcard $(ANA)/SAXS_12id/reduced_data/*.csv)
RAD_DATA  := $(wildcard $(ANA)/Rad_Dam_Check/cluster_results/*.hdf)
COMMON    := $(ANA)/common/acs_style.py

.PHONY: all figures papers main si cover check clean distclean env

all: figures papers

# --- environment ------------------------------------------------------------
env:
	conda env create -f $(ANA)/environment.yml

# --- figures ----------------------------------------------------------------
figures: $(FIGURES)
	@echo "all figures up to date"

$(F_12ID): $(ANA)/SAXS_12id/Plot_12ID.py $(CSV_DATA) $(COMMON)
	cd $(ANA)/SAXS_12id && $(PY) Plot_12ID.py

$(F_GUINIER): $(ANA)/SAXS_12id/Guinier_Plot.py $(CSV_DATA) $(COMMON)
	cd $(ANA)/SAXS_12id && $(PY) Guinier_Plot.py

# one script, three figures
$(F_XPCS) &: $(ANA)/SAXPCS_8id/saxpcs.py $(ANA)/SAXPCS_8id/abs_xsec.py $(XPCS_DATA) $(COMMON)
	cd $(ANA)/SAXPCS_8id && $(PY) saxpcs.py

$(F_EVOL): $(ANA)/SAXPCS_8id/saxs_evolution.py $(ANA)/SAXPCS_8id/abs_xsec.py $(XPCS_DATA) $(COMMON)
	cd $(ANA)/SAXPCS_8id && $(PY) saxs_evolution.py

$(F_GRID): $(ANA)/SAXPCS_8id/g2_grid_SI.py $(XPCS_DATA) $(COMMON)
	cd $(ANA)/SAXPCS_8id && $(PY) g2_grid_SI.py

$(F_RAD): $(ANA)/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py $(ANA)/Rad_Dam_Check/Flux_Cal.py $(RAD_DATA) $(COMMON)
	cd $(ANA)/Rad_Dam_Check && $(PY) g2_SAXPCS_Rad_Cali.py

# --- documents --------------------------------------------------------------
# bibtex runs only when the .aux actually declares a bibliography, so the
# cover letter (which has none) does not fail the build.
define compile
	@mkdir -p $(BUILD) $(OUT)
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@cd $(MS) && if grep -q '\\bibdata' build/$(1).aux 2>/dev/null; then \
	    bibtex build/$(1) >/dev/null; \
	    pdflatex -interaction=nonstopmode -halt-on-error \
	        -file-line-error -output-directory=build $(1).tex >/dev/null; \
	fi
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@cp $(BUILD)/$(1).pdf $(OUT)/$(2)
	@echo "  -> manuscript/examples/$(2) ($$(pdfinfo $(OUT)/$(2) | awk '/^Pages/{print $$2}') pages)"
endef

papers: main si cover

main: $(FIGURES) $(MS)/main.tex $(MS)/reference.bib
	$(call compile,main,main.pdf)

si: $(FIGURES) $(MS)/si.tex $(MS)/reference.bib
	$(call compile,si,SI.pdf)

cover: $(MS)/cover_letter.tex $(MS)/cover_letter_header.tex
	$(call compile,cover_letter,cover_letter.pdf)

# --- verification -----------------------------------------------------------
# NB: LaTeX writes Citation `key' undefined -- backtick then straight quote.
check:
	@echo "== undefined citations / references =="
	@grep -hoE "(Citation|Reference) \`[^']+' undefined" \
	    $(BUILD)/main.log $(BUILD)/si.log 2>/dev/null | sort -u || true
	@echo "== bibtex database/entry problems =="
	@grep -hE "couldn't open|I didn't find a database entry" \
	    $(BUILD)/*.blg 2>/dev/null | sort -u | head -20 || true
	@echo "== missing figures =="
	@grep -hE "No file|not found" $(BUILD)/*.log 2>/dev/null | grep -iE '\.pdf|\.png' || true
	@echo "== done (no lines under a heading means clean) =="

# --- cleaning ---------------------------------------------------------------
clean:
	rm -f $(BUILD)/*.aux $(BUILD)/*.log $(BUILD)/*.out $(BUILD)/*.bbl \
	      $(BUILD)/*.blg $(BUILD)/*.toc $(BUILD)/*.spl $(BUILD)/*.fls

distclean:
	rm -rf $(BUILD)
