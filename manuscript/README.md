# Manuscript sources

LaTeX sources for the Nano Letters submission. Build everything from the
repository root:

```
make papers     # -> build/main.pdf, build/si.pdf, build/cover_letter.pdf
make check      # report undefined citations/references or missing figures
```

| Source | Class | Output |
|---|---|---|
| `main.tex` | `achemso` (`journal=nalefd`, Nano Letters) | `build/main.pdf` |
| `si.tex` | `article` + `natbib`, `achemso` bib style | `build/si.pdf` |
| `cover_letter.tex` | `article` + `cover_letter_header.tex` | `build/cover_letter.pdf` |

Intermediate files go to `build/` (gitignored); only the finished PDFs are
written to `build/`, together with the LaTeX aux files. Source -- the `.tex`
files, `reference.bib`, `figures/` and the vendored `texmf/` tree -- stays
outside `build/`, because `make distclean` deletes that directory wholesale.

## Where the figures come from

`\graphicspath{{figures/}}` is left as-is and the `.tex` files use short figure
names. The build puts the data-derived figure directories on `TEXINPUTS`, so
those names resolve even though the files live under `analysis/`:

| Figure name in `.tex` | Actual location |
|---|---|
| `Setup.pdf` (Fig. 1) | `manuscript/figures/` |
| `FigureS1_Sample_Cells.png`, `VideoS1_Still.png` | `manuscript/figures/` |
| `Figure2_SAXS_WAXS.pdf`, `FigureS5_Guinier.pdf` | `analysis/SAXS_12id/` |
| `Figure3_Isothermal_SAXPCS.pdf`, `FigureS3_Calibration.pdf`, `FigureS6_SAXS_Evolution.pdf`, `FigureS7_Fit_Parameters.pdf`, `FigureS8_g2_Grid.pdf` | `analysis/SAXPCS_8id/` |
| `FigureS4_Flux_Control.pdf` | `analysis/Rad_Dam_Check/` |

Figure 1 is the approved illustrator schematic; everything else with a `.pdf`
extension is regenerated from the experimental data by `make figures`.

## texmf/

`achemso` (the ACS class, not in the system TeX Live here) and `captdef` (needed
by `cover_letter_header.tex`) are vendored into `texmf/`, which the Makefile
exposes as `TEXMFHOME`. This keeps the build self-contained — no root install
and no `tlmgr` step. `texmf/tex/latex/achemso/LICENSE.md` carries the upstream
LGPL/LPPL terms.

Two BibTeX details worth knowing if you build by hand rather than via `make`:

- `achemso` writes a control file `acs-main.bib` next to the `.aux`. Because the
  build uses `-output-directory=build`, that directory must be on `BIBINPUTS`
  or BibTeX silently drops *every* database and all citations come out
  undefined.
- `main.tex` and `si.tex` both need three `pdflatex` passes around one `bibtex`
  pass to settle references.
