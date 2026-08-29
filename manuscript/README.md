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

Everything LaTeX produces goes to `build/`: the finished PDFs, which are kept
in git so collaborators can read them without a TeX installation, and the aux
files, which `manuscript/.gitignore` excludes (`build/*` with `!build/*.pdf`). Source -- the `.tex` files,
`reference.bib`, `figures/` and the vendored `texmf/` tree -- stays outside
`build/`, because `make distclean` deletes that directory wholesale.

## Where the figures come from

`\graphicspath{{figures/}}` is left as-is and the `.tex` files use short figure
names. The build puts the data-derived figure directories on `TEXINPUTS`, so
those names resolve even though the files live under `analysis/`:

| Figure name in `.tex` | Actual location |
|---|---|
| `Setup.pdf` (Fig. 1) | `manuscript/figures/` |
| `FigureS2_Sample_Cells.png` (Fig. S2), `VideoS1_Still.png` | `manuscript/figures/` |
| `Figure2_SAXS_WAXS.pdf`, `FigureS7_Guinier.pdf` | `analysis/SAXS_12id/` |
| `Figure3_Isothermal_SAXPCS.pdf`, `FigureS3_Calibration.pdf`, `FigureS4_Contrast.pdf`, `FigureS6_Thermal_Cycle.pdf`, `FigureS8_SAXS_Evolution.pdf`, `FigureS9_Fit_Parameters.pdf`, `FigureS10_g2_Grid.pdf` | `analysis/SAXPCS_8id/` |
| `FigureS5_Flux_Control.pdf` | `analysis/Rad_Dam_Check/` |

Figure S1 is a LaTeX-typeset sequence box, so it has no image file; every other
figure file is named for the number it carries in the text.

Figure 1 is the approved illustrator schematic; everything else with a `.pdf`
extension is regenerated from the experimental data by `make figures`.

## Referring to figures

No figure number is ever typed by hand. Each figure carries a label --
`fig:setup`, `fig:static`, `fig:xpcs` in `main.tex`, `fig:S1`--`fig:S10` in
`si.tex` -- and the text writes `Figure~\ref{fig:S6}`, appending any panel
letter literally (`Figure~\ref{fig:S6}b`). The two documents also cite each
other's figures, which `\ref` alone cannot reach, so each loads `xr` and
imports the other under a prefix: `\externaldocument[SI-]{build/si}` in
`main.tex` and `\externaldocument[MT-]{build/main}` in `si.tex`. A reference
across documents therefore reads `Figure~\ref{SI-fig:S7}` or
`Figure~\ref{MT-fig:xpcs}a`. Because that makes each document depend on the
other's `.aux`, `make papers` runs si -> main -> si; `make check` reports
anything left unresolved, which is what a stale or missing `.aux` looks like.
Renumbering or reordering figures then needs no edit outside the labels.

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
