# Build everything for the ELP gelation paper: figures from data, then PDFs.
#
#   make env        create the conda environment (once)
#   make figures    regenerate every data-derived figure in analysis/
#   make papers     compile main.pdf, si.pdf, cover_letter.pdf into manuscript/build/
#   make all        figures + papers
#   make check      report undefined citations / references / missing figures
#   make robustness recompute the SI 5.3 robustness numbers (no figure)
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

CONDA_ENV := env-2021-LLDP
PY        := conda run --no-capture-output -n $(CONDA_ENV) python

# achemso (ACS/Nano Letters class) and captdef are not in the system TeX Live,
# so they live in a repo-local texmf tree -- this keeps the build self-contained.
export TEXMFHOME := $(MS)/texmf
# Bare figure names in the .tex resolve through here; that is what lets the
# manuscript keep short names while the files live under analysis/.
export TEXINPUTS := .:$(MS)/figures:$(ANA)/SAXPCS_8id:$(ANA)/Rad_Dam_Check:
# $(BUILD) must be on BIBINPUTS: achemso writes a control file (acs-main.bib)
# next to the .aux, and BibTeX aborts all databases if it cannot open it.
export BIBINPUTS := .:$(MS):$(BUILD):
export MPLBACKEND := Agg

# --- data-derived figures: script -> outputs --------------------------------
F_RAD     := $(ANA)/Rad_Dam_Check/FigureS6_Flux_Control.pdf
F_TSERIES := $(ANA)/SAXPCS_8id/Figure3_Isothermal_SAXPCS.pdf \
             $(ANA)/SAXPCS_8id/FigureS7_Temperature_2800.pdf \
             $(ANA)/SAXPCS_8id/FigureS8_Temperature_2825.pdf \
             $(ANA)/SAXPCS_8id/FigureS9_Temperature_2850.pdf \
             $(ANA)/SAXPCS_8id/FigureS10_Temperature_2875.pdf
F_ARREST  := $(ANA)/SAXPCS_8id/Figure4_Arrest_vs_Temperature.pdf \
             $(ANA)/SAXPCS_8id/FigureS11_Collapse_AllQ.pdf

# The 2022-1 arrest analysis was retired with the move to the 2021-3
# temperature series: saxs_evolution.py, g2_grid_SI.py and the fit-parameter
# panels of saxpcs.py no longer have a figure in either document.  saxpcs.py
# is still built for the absolute-scale calibration, which both beamtimes share.
FIGURES := $(F_RAD) $(F_TSERIES) $(F_ARREST)

RAD_DATA  := $(wildcard $(ANA)/Rad_Dam_Check/cluster_results/*.hdf)
TSERIES_DATA := $(wildcard $(ANA)/SAXPCS_8id/data_2021-3/*.hdf) \
                $(ANA)/SAXPCS_8id/data_2021-3/trace_2021-3.csv
COMMON    := $(ANA)/common/acs_style.py
XFIT      := $(ANA)/SAXPCS_8id/xpcs_fit.py
# nexus_read.py is imported by the four 8-ID-I figure scripts below, so a change
# to it has to rebuild all four; without it here they would go quietly stale.
NEXUS     := $(ANA)/SAXPCS_8id/nexus_read.py

.PHONY: all figures papers main si cover check robustness clean distclean env

all: figures papers

# --- environment ------------------------------------------------------------
env:
	conda env create -f $(ANA)/environment.yml

# --- figures ----------------------------------------------------------------
figures: $(FIGURES)
	@echo "all figures up to date"

# one script, three figures

$(F_TSERIES) &: $(ANA)/SAXPCS_8id/temperature_series_2021-3.py $(ANA)/SAXPCS_8id/abs_xsec.py \
                $(TSERIES_DATA) $(COMMON) $(XFIT) $(NEXUS)
	cd $(ANA)/SAXPCS_8id && $(PY) temperature_series_2021-3.py

# the four 2021-3 holds
$(F_ARREST) &: $(ANA)/SAXPCS_8id/arrest_vs_temperature.py \
             $(TSERIES_DATA) $(COMMON) $(XFIT) $(NEXUS)
	cd $(ANA)/SAXPCS_8id && $(PY) arrest_vs_temperature.py

$(F_RAD): $(ANA)/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py $(ANA)/Rad_Dam_Check/Flux_Cal.py $(RAD_DATA) $(COMMON)
	cd $(ANA)/Rad_Dam_Check && $(PY) g2_SAXPCS_Rad_Cali.py

# --- documents --------------------------------------------------------------
# bibtex runs only when the .aux actually declares a bibliography, so the
# cover letter (which has none) does not fail the build.
# Everything LaTeX produces -- aux files and the finished PDFs alike -- lands in
# build/.  Source (.tex, .bib, figures/, texmf/) stays outside it, because
# `make distclean` deletes build/ wholesale.
define compile
	@mkdir -p $(BUILD)
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@cd $(MS) && if grep -q '\\bibdata' build/$(1).aux 2>/dev/null; then \
	    bibtex build/$(1) >/dev/null; \
	    pdflatex -interaction=nonstopmode -halt-on-error \
	        -file-line-error -output-directory=build $(1).tex >/dev/null; \
	fi
	@cd $(MS) && pdflatex -interaction=nonstopmode -halt-on-error \
	    -file-line-error -output-directory=build $(1).tex >/dev/null
	@echo "  -> manuscript/build/$(1).pdf ($$(pdfinfo $(BUILD)/$(1).pdf | awk '/^Pages/{print $$2}') pages)"
endef

# The xr package makes each document read the other's .aux for cross-document figure
# numbers, so on a clean tree the pair has to be run si -> main -> si before
# both directions resolve.  `make main` / `make si` alone stay available for
# iterating on one document; `make check` reports anything still unresolved.
papers: $(FIGURES) $(MS)/main.tex $(MS)/si.tex $(MS)/reference.bib \
        $(MS)/cover_letter.tex $(MS)/cover_letter_header.tex
	$(call compile,si)
	$(call compile,main)
	$(call compile,si)
	$(call compile,cover_letter)

main: $(FIGURES) $(MS)/main.tex $(MS)/reference.bib
	$(call compile,main)

si: $(FIGURES) $(MS)/si.tex $(MS)/reference.bib
	$(call compile,si)

cover: $(MS)/cover_letter.tex $(MS)/cover_letter_header.tex
	$(call compile,cover_letter)

# --- verification -----------------------------------------------------------
# NB: LaTeX writes  Reference `key' on page N undefined on input line M.
# The "undefined" does NOT follow the closing quote directly, so the pattern
# has to allow the page clause in between.  An earlier version anchored on
# "' undefined" and therefore matched nothing, reporting every build clean.
check:
	@echo "== undefined citations / references =="
	@grep -hoE "(Citation|Reference) \`[^']+'.*undefined" \
	    $(BUILD)/main.log $(BUILD)/si.log 2>/dev/null | sort -u || true
	@echo "== bibtex database/entry problems =="
	@grep -hE "couldn't open|I didn't find a database entry" \
	    $(BUILD)/*.blg 2>/dev/null | sort -u | head -20 || true
	@echo "== missing figures =="
	@grep -hE "No file|not found" $(BUILD)/*.log 2>/dev/null | grep -iE '\.pdf|\.png' || true
	@echo "== done (no lines under a heading means clean) =="

# --- the SI 5.3 robustness check --------------------------------------------
# Numbers only; this one draws nothing, so it is a prerequisite of nothing and
# `make all` does not run it.  It refits all 34 groups three ways over a p_fast
# grid and takes roughly half an hour.  Run it when SI Section 5.3 is rewritten
# or the 2021-3 data are reprocessed, and paste its output into that section.
robustness:
	cd $(ANA)/SAXPCS_8id && $(PY) general_model_check.py

# --- cleaning ---------------------------------------------------------------
clean:
	rm -f $(BUILD)/*.aux $(BUILD)/*.log $(BUILD)/*.out $(BUILD)/*.bbl \
	      $(BUILD)/*.blg $(BUILD)/*.toc $(BUILD)/*.spl $(BUILD)/*.fls

# NB: removes the compiled PDFs too.  They are regenerable by `make papers`.
distclean:
	rm -rf $(BUILD)
