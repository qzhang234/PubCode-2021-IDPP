# Manuscript sources

LaTeX sources for the PNAS submission. Build everything from the repository
root:

```
make papers     # -> build/main.pdf, build/si.pdf, build/cover_letter.pdf
make check      # report undefined citations/references or missing figures
```

| Source | Class | Output |
|---|---|---|
| `main.tex` | `achemso` (PNAS accepts any format at initial submission) | `build/main.pdf` |
| `si.tex` | `article` + `natbib`, `achemso` bib style | `build/si.pdf`, the SI Appendix |
| `cover_letter.tex` | `article` + `cover_letter_header.tex` | `build/cover_letter.pdf` |

Everything LaTeX produces goes to `build/`: the finished PDFs, which are kept
in git so collaborators can read them without a TeX installation, and the aux
files, which `manuscript/.gitignore` excludes. Source -- the `.tex` files,
`reference.bib`, `figures/` and the vendored `texmf/` tree -- stays outside
`build/`, because `make distclean` deletes that directory wholesale.

## Where the figures come from

The `.tex` files use short figure names, and the build puts the data-derived
figure directories on `TEXINPUTS`, so those names resolve even though the files
live under `analysis/`:

| Figure name in `.tex` | Actual location |
|---|---|
| `Setup.pdf` (Fig. 1), `FigureS1b_SDS_PAGE.png` | `manuscript/figures/` |
| `Figure3_Isothermal_SAXPCS.pdf`, `Figure4_Arrest_vs_Temperature.pdf`, `FigureS2_Calibration.pdf`, `FigureS7`-`S10_Temperature_*.pdf`, `FigureS11_Collapse_AllQ.pdf`, `FigureS12_Fast_Mode.pdf` | `analysis/SAXPCS_8id/` |
| `FigureS6_Flux_Control.pdf` | `analysis/Rad_Dam_Check/` |

The file names carry an earlier numbering; the printed numbers come from the
labels. Figure 1 is the approved illustrator schematic; every other `.pdf`
figure is regenerated from the experimental data by `make figures`.

## Referring to figures

No figure number is ever typed by hand. Each figure carries a label, and the
text cites it through `\ref`, appending any panel letter literally. `main.tex`
uses the PNAS forms, `Fig.~\ref{fig:xpcs}c` and `Eq.~\ref{eq:twomode}`;
`si.tex` keeps `Figure~\ref{...}`. The two documents cite each other's figures,
so each loads `xr` and imports the other under a prefix:
`\externaldocument[SI-]{build/si}` in `main.tex` and
`\externaldocument[MT-]{build/main}` in `si.tex`. A reference across documents
therefore reads `SI Appendix, Fig.~\ref{SI-fig:collapse}` in the main text, or
`Figure~\ref{MT-fig:xpcs}a` in the SI.
Because each document then depends on the other's `.aux`, `make papers` runs
si -> main -> si; `make check` reports anything left unresolved.

## texmf/

`achemso` (not in the system TeX Live here) and `captdef` (needed by
`cover_letter_header.tex`) are vendored into `texmf/`, which the Makefile
exposes as `TEXMFHOME`, so the build needs no root install.
`texmf/tex/latex/achemso/LICENSE.md` carries the upstream LGPL/LPPL terms.

Two BibTeX details if you build by hand rather than via `make`:

- `achemso` writes a control file `acs-main.bib` next to the `.aux`. Because
  the build uses `-output-directory=build`, that directory must be on
  `BIBINPUTS` or BibTeX silently drops every database.
- `main.tex` and `si.tex` both need three `pdflatex` passes around one `bibtex`
  pass to settle references.
