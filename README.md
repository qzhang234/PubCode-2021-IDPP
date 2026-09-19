# Reversible β-Sheet-like Ordering and Dynamic Arrest in an Elastin-like Polypeptide

Everything needed to reproduce the Nano Letters submission: the reduced
experimental data, the analysis code that turns it into every figure, and the
LaTeX sources for the manuscript, Supporting Information and cover letter.

The measurements combine SAXS/WAXS at APS beamline **12-ID-B** (static
structure across the thermal transition) with SA-XPCS at **8-ID-I** (dynamics
during an isothermal hold at 30 °C).

```
make env        # create the conda environment, once
make figures    # regenerate every data-derived figure from the reduced data
make papers     # compile main.pdf, si.pdf, cover_letter.pdf into manuscript/build/
make all        # figures + papers
make check      # report undefined citations / references / missing figures
make clean      # remove LaTeX aux files, keep the PDFs
```

## Reproducibility

Two scripts, and only two, read raw beamline storage:

| Script | Reads | Writes |
|---|---|---|
| `analysis/SAXS_12id/Read_12ID_SAWAXS.py` | 12-ID-B `.avg` profiles | `analysis/SAXS_12id/reduced_data/Merged_*.csv` |
| `analysis/SAXPCS_8id/average_ranges.py` | 8-ID-I reprocessed NeXus results | `analysis/SAXPCS_8id/data/` |

Both of their outputs are committed, so **`make figures` and `make papers` run
from a clone of this repository alone** — no beamline account, no network, no
raw data. The two scripts above have to be rerun only if the raw data are
reprocessed, and they will only run at the beamline.

### Code that reads outside this repository

There are exactly three places in the Python where a path can point outside the
repository, and **none of them is on the path from the committed data to a
figure** — `make figures` never executes any of them.

| Where | Out-of-repo location | How the path is set |
|---|---|---|
| `analysis/SAXS_12id/Read_12ID_SAWAXS.py`, `fn_path` (line 54) | `/home/8-id-i/2021-1/12-id-b/ZuoApr13/Processed/` — the 12-ID-B reduced `.avg` profiles from the 2021-1 run | hard-coded |
| `analysis/SAXPCS_8id/average_ranges.py`, `prefix` (line 80) | `/home/8-id-i/2022-1/babnigg202203_nexus/reprocess_results` — the 8-ID-I per-acquisition cluster results from the 2022-1 run | hard-coded |
| `analysis/common/utils.py`, beamline-side helpers | wherever the caller points them | **no default** — `prefix` must be passed in |

The `utils.py` helpers in question are `read_keys_from_files`,
`read_keys_from_files_parallel`, `get_temperature`,
`read_temperature_from_files`, `_read_xpcs_hdf`, `process_group` and
`process_group_by_range`. They are kept for provenance and for reprocessing at
the beamline; two of them also import `pyxpcsviewer`, which is deliberately not
in `environment.yml`. The only two functions in that module that a figure
depends on are `outlier_removal` and `average_datasets`, and both take numpy
arrays — neither ever opens a file.

Those three are the only absolute path literals in the whole tree. Every other
script anchors its input and output on
`os.path.dirname(os.path.abspath(__file__))`, so it reads and writes inside its
own directory and runs correctly from any working directory.

## Source data integrity

No data file is modified by the act of producing a figure. Checked three ways.

**1. Every read is read-only.** Every `h5py.File` call in every plotting script
opens with mode `'r'`. The tree contains exactly one `'r+'` open, in
`average_ranges.py:save_average`, and it acts on the *destination* of a
`shutil.copyfile` under `analysis/SAXPCS_8id/data/` — never on a source file.
The complete list of things any script writes is: the four
`SAXS_12id/reduced_data/*.csv`, the contents of `SAXPCS_8id/data/`, and the ten
figure PDFs. One in-place-mutation hazard is guarded explicitly:
`apply_cross_corr_threshold` copies its input before rewriting NaNs and
non-positive values to make the similarity metric well defined, so those
injected values cannot leak back into `saxs_1d`.

**2. Syscall trace.** A full `make figures` was run under `strace`, tracing
`open`/`openat`/`unlink`/`rename`/`truncate` across all child processes. The
complete set of paths opened for writing is the ten figure PDFs, four
`__pycache__/*.pyc` files, and matplotlib's scratch files in `/tmp`. No `.hdf`
or `.csv` is opened for writing, and no `unlink`, `rename` or `truncate` touches
a data file. All 46 data files appear in the trace as `O_RDONLY` only. Nothing
outside the repository is read except Python/conda libraries, matplotlib's font
cache and `/proc`–`/sys`.

**3. Checksums.** md5 of all 46 committed data files — 40 `.hdf` and 6 `.csv` —
taken before and after a full `make figures`, repeated. All 46 identical every
time.

Related, and for the same reason: nothing is rescaled between the calibration
and the plot. `abs_xsec_coef()` is built from millimetres, so every absolute
intensity axis in this work is **mm⁻¹**, plotted as the coefficient comes out
with no conversion factor applied. (The literature more often quotes cm⁻¹, which
is ten times larger.) The single non-unity coefficient anywhere in the analysis
is the empirical `alpha_WA = 0.95` in the 12-ID WAXS reduction, disclosed in SI
Section 3.1.

## Layout

| Path | What it is |
|---|---|
| `Makefile` | The whole build: figures, the three PDFs, and `make check`. |
| `README.md` | This file. |
| `ACS_author_guide.pdf` | ACS/Nano Letters author guidelines. Appendix 2 is the source of every number in `analysis/common/acs_style.py`; the Editorial Policies section is the source of the length caps. |
| `analysis/` | Data reduction and figure generation. See `analysis/README.md`. |
| `manuscript/` | LaTeX sources and compiled PDFs. See `manuscript/README.md`. |
| `.gitignore` | Excludes Claude Code session state (`.claude/`, `CLAUDE.md`, `claude_log/`). |
| `.loglogin` | Beamline login log, carried over from the acquisition account. Not used by anything here. |

### `analysis/`

| Path | What it is |
|---|---|
| `environment.yml` | Conda environment `env-2021-LLDP`: python 3.11, numpy, scipy, matplotlib 3.9.4 (pinned for byte-stable figures), h5py, pandas. |
| `common/acs_style.py` | Single source of ACS figure compliance: 3.33 / 7.0 in column widths, 8 pt Arial everywhere, ≥ 0.5 pt lines, TrueType embedding, `save_fig()`. Imported by every plotting script. |
| `common/utils.py` | `outlier_removal` and `average_datasets`, used by `average_ranges.py`; plus beamline-side reduction helpers kept for provenance that cannot run from a clone. |
| `SAXS_12id/` | 12-ID-B: `Read_12ID_SAWAXS.py` (reduction), `Plot_12ID.py` → **Fig. 2**, `Guinier_Plot.py` → **Fig. S7**, and the committed `reduced_data/*.csv`. |
| `SAXPCS_8id/` | 8-ID-I: `average_ranges.py` (reduction) and the committed `data/`; `saxpcs.py` → **Figs. 3, S3, S9**; `contrast_calibration.py` → **Fig. S4**; `thermal_cycle.py` → **Fig. S6**; `saxs_evolution.py` → **Fig. S8**; `g2_grid_SI.py` → **Fig. S10**. `abs_xsec.py` holds the absolute-cross-section calibration and `xpcs_fit.py` the shared two-mode g2 model, so every figure uses one implementation of each. `timelist_2022-1.txt` recovers acquisition times the 2025 reprocessing overwrote; `nexus_manual.txt` documents the NeXus layout. |
| `Rad_Dam_Check/` | Flux-dependence control → **Fig. S5**, with its four result files, the PIN-diode calibration sheet and certificate. See `analysis/Rad_Dam_Check/README.md`. |

### `manuscript/`

| Path | What it is |
|---|---|
| `main.tex` | The Letter (`achemso`, `journal=nalefd`). |
| `si.tex` | Supporting Information, a separate document as Nano Letters requires. |
| `cover_letter.tex`, `cover_letter_header.tex` | Cover letter and its preamble. |
| `reference.bib` | Shared bibliography for both documents. |
| `figures/Setup.pdf` | **Fig. 1**, the approved schematic. Not generated by any script; do not regenerate. |
| `figures/FigureS1_Sample_Cells.png` | **Fig. S2** (the file name predates the current SI figure order). Fig. S1 is a LaTeX-typeset sequence box with no image file. |
| `figures/VideoS1_Still.png`, `figures/anl.png` | Video S1 still frame and the Argonne logo. |
| `texmf/` | Vendored `achemso` and `captdef`, exposed as `TEXMFHOME` so the build needs no root TeX install. |
| `build/` | Build output. The three finished PDFs are committed; the LaTeX aux files are ignored. |

Figures are named exactly as the `.tex` files reference them and are embedded at
natural size, so 8 pt drawn is 8 pt printed. The build resolves the names
through `TEXINPUTS` rather than copying the files, so there is one copy of each
figure and it is the one the script wrote.

## Submission targets

Nano Letters Editorial Policies: under 3000 words, no more than 5 figures,
abstract no more than 150 words, no subsection titles, Supporting Information as
a separate file.
