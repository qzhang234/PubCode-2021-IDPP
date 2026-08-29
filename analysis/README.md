# llps-saxpcs

Data reduction and figure-generation code for the manuscript
"Reversible β-Sheet-like Ordering and XPCS-Resolved Dynamic Arrest during
Elastin-like Polypeptide Phase Separation" (main text + Supplemental Information).

The experiment combines small/wide-angle X-ray scattering (SAXS/WAXS) at
12-ID-B with X-ray photon correlation spectroscopy (XPCS) at 8-ID-I to follow
an elastin-like polypeptide (ELP) sample through a liquid-liquid phase
separation (LLPS): SAXS/WAXS resolves the static structure (Guinier radius,
the emergent β-sheet-like correlation peaks), while XPCS resolves the
dynamics (g2 relaxation) as the sample ages after a temperature jump. This
repo goes from raw beamline files to the PDF figures used in the paper.

## Repository layout

- **`common/`** — code shared by more than one plotting script.
  - `utils.py`: cross-correlation-based outlier removal and averaging for XPCS
    cluster-result files. Only `outlier_removal` and `average_datasets` are on
    the path to a figure (both used by `average_ranges.py`); the rest of the
    module is beamline-side tooling that walks raw result directories and
    cannot run from a clone of this repository.
  - `acs_style.py`: the shared ACS figure-formatting module (column widths,
    fonts, line/marker sizes, panel-letter helper) used by every script below
    — see [Figure formatting](#figure-formatting).
- **`SAXS_12id/`** — 12-ID-B SAXS/WAXS: Guinier analysis and the merged
  SAXS+WAXS structure-factor plot.
- **`SAXPCS_8id/`** — 8-ID-I XPCS: g2 relaxation, fit parameters, absolute
  SAXS cross-section evolution, and the beamline flux/photon calibration.
- **`Rad_Dam_Check/`** — radiation-damage control: confirms the XPCS dynamics
  are flux-independent over the attenuation range used in the main
  measurement.

## Data flow

### `SAXS_12id/`

```
raw .avg files (12-ID-B beamline storage; BEAMLINE-ONLY)
    --> Read_12ID_SAWAXS.py   (merge SAXS+WAXS, scale WAXS onto SAXS, propagate errors)
    --> reduced_data/Merged_*.csv
    --> Guinier_Plot.py       --> FigureS7_Guinier.pdf
    --> Plot_12ID.py          --> Figure2_SAXS_WAXS.pdf
```

`Read_12ID_SAWAXS.py` is one of only two scripts in the repository that read
raw beamline storage (the other is `SAXPCS_8id/average_ranges.py`). Its outputs
are committed under `reduced_data/`, so it only has to be rerun if the raw
12-ID-B data are reprocessed.

`Guinier_Plot.py` fits `ln[I(Q)]` vs. `Q^2` in the low-Q Guinier region of the
10 °C reference to get R_g. `Plot_12ID.py` overlays the merged reference and
measurement curves, annotates the high-Q power law, and fits the two
β-sheet-like correlation peaks with Gaussians in an inset.

### `SAXPCS_8id/`

```
raw cluster-result .hdf files (8-ID-I, per-file g2, g2_err, saxs_1d, IC readout)
    --> average_ranges.py    --> data/
data/
    --> saxpcs.py              --> Figure3_Isothermal_SAXPCS.pdf,
                                    FigureS9_Fit_Parameters.pdf,
                                    FigureS3_Calibration.pdf
    --> saxs_evolution.py      --> FigureS8_SAXS_Evolution.pdf
    --> g2_grid_SI.py          --> FigureS10_g2_Grid.pdf
    --> thermal_cycle.py       --> FigureS6_Thermal_Cycle.pdf
    --> contrast_calibration.py --> FigureS4_Contrast.pdf
```

`average_ranges.py` is the second of the two beamline-only scripts and the
only one in this directory that reads the beamline storage. It reduces
them into `data/`, which is committed: an averaged NeXus file per frame range
(outlier-filtered, carrying its own mean ion-chamber readings and the list of
acquisitions that survived), the per-acquisition temperature trace of the
thermal-cycling sequence, and the 50 individual correlation functions of the
contrast standard. Every figure script below reads only `data/`, so the whole
analysis can be rerun from this repository alone.

One group is not reduced purely automatically. `MANUAL_EXCLUDE` at the top of
`average_ranges.py` lists acquisitions dropped by hand on top of
`outlier_removal()`, which cuts on the *shape* of log10 I(q) and therefore
cannot see a curve that spikes in the three lowest q bins and is normal
everywhere else. The table currently holds eleven acquisitions of B0083, cycle
4 of the Figure S6 series; SI Section 7.1 states what removing them does and
what a uniform cut across all seven cycles would give instead. Re-reduce just
that group with `python average_ranges.py B0083`.

`saxpcs.py` reads the averaged B0147 (isothermal 30 °C) and D0138 (buffer)
HDF files directly with h5py, fits a double-exponential (KWW) model to g2(τ)
globally over the five lowest q bins at once, with the two stretching exponents
shared across q, then fits each elapsed time's τ_fast(Q) and τ_slow(Q)
to a power law τ = A·Q^γ (weighted log-log regression with error propagation)
to extract the scaling exponents γ_fast and γ_slow. It also converts SAXS
intensities to an absolute differential cross section via `abs_xsec.py`
(per-file ion-chamber-based coefficients — see the calibration derivation in
that file), and calibrates the upstream ion-chamber reading against the
incident photon flux. `saxs_evolution.py`
reuses the same absolute-cross-section calibration to show the full SAXS
evolution across all B0147 files. `g2_grid_SI.py` repeats the g2 fit for the
four secondary q bins not shown in the main figure. `thermal_cycle.py` puts the
seven thermal cycles of the reversibility control on the same absolute scale
with the same buffer subtraction, and `contrast_calibration.py` measures the
instrumental speckle contrast β on the static glass standard that every g2 fit
then holds fixed.

### `Rad_Dam_Check/`

```
cluster-result .hdf files at 4 attenuation settings (E0171-E0174)
Pind_calibration.csv (PIN-diode <-> ion-chamber calibration sheet)
    --> g2_SAXPCS_Rad_Cali.py --> FigureS5_Flux_Control.pdf
```

Converts each measurement's ion-chamber reading to an on-sample photon flux
(the nominal attenuator labels are not monotonic in flux, so the true flux is
recovered from the calibration sheet) and shows that both the SAXS I(Q) and
the XPCS g2(τ) are unchanged across the flux range used in the main
measurement — i.e., the reported dynamics are not a radiation-damage
artifact. `Flux_Cal.py` holds the underlying flux-conversion formula. The four
result files are committed under `cluster_results/`, so this figure also
rebuilds from the repository alone.

## Figure-to-manuscript mapping

Each output PDF is named exactly as `main.tex` / `si.tex` reference it, so the
LaTeX build picks the files up straight from this directory — only the path
differs, and the top-level `Makefile` resolves it via `TEXINPUTS`.

| Output PDF | Script | Figure | Panels |
|---|---|---|---|
| `SAXS_12id/Figure2_SAXS_WAXS.pdf` | `Plot_12ID.py` | Main Fig. 2 | 1 (+ inset) |
| `SAXS_12id/FigureS7_Guinier.pdf` | `Guinier_Plot.py` | SI Fig. S7 | 1 |
| `SAXPCS_8id/Figure3_Isothermal_SAXPCS.pdf` | `saxpcs.py` | Main Fig. 3 | 3 |
| `SAXPCS_8id/FigureS9_Fit_Parameters.pdf` | `saxpcs.py` | SI Fig. S9 | 4 (2×2) |
| `SAXPCS_8id/FigureS3_Calibration.pdf` | `saxpcs.py` | SI Fig. S3 | 2 |
| `SAXPCS_8id/FigureS4_Contrast.pdf` | `contrast_calibration.py` | SI Fig. S4 | 2 |
| `SAXPCS_8id/FigureS6_Thermal_Cycle.pdf` | `thermal_cycle.py` | SI Fig. S6 | 3 |
| `SAXPCS_8id/FigureS8_SAXS_Evolution.pdf` | `saxs_evolution.py` | SI Fig. S8 | 1 |
| `SAXPCS_8id/FigureS10_g2_Grid.pdf` | `g2_grid_SI.py` | SI Fig. S10 | 4 (2×2) |
| `Rad_Dam_Check/FigureS5_Flux_Control.pdf` | `g2_SAXPCS_Rad_Cali.py` | SI Fig. S5 | 2 |

Main Fig. 1 (`manuscript/figures/Setup.pdf`) is an illustrator schematic and is
not produced by any script here. SI Fig. S1 is a LaTeX-typeset box of
oligonucleotide sequences with no image file. SI Fig. S2 is a photograph of the
sample cells; its image file is `manuscript/figures/FigureS1_Sample_Cells.png`,
whose name predates the current SI figure order. The Video S1 still frame
(`manuscript/figures/VideoS1_Still.png`) is likewise a photograph.

## Environment

`environment.yml` pins the five scientific packages the figure scripts need:

```
conda env create -f environment.yml     # creates env-2021-LLDP
conda activate env-2021-LLDP
```

matplotlib is pinned to 3.9.4 for reproducible layout: `tight_layout` spacing
metrics shift slightly between minor releases, which moves axes, keys, and
annotations inside the frame. The page size itself no longer moves with the
version — `save_fig()` writes at the exact figure size with no
`bbox_inches='tight'` trimming — but the pin preserves the collision-free
placements checked for this submission. `pyxpcsviewer` is *not* required — it is imported
lazily inside `common/utils.py` helpers that read raw beamline files, which is
not on the path from the reduced data in this repo to any figure.

Regenerate everything from the repository root with `make figures`, or run a
single script from its own directory (`cd SAXPCS_8id && python saxpcs.py`).

## Figure formatting

All figures follow the ACS "Preparing Graphics" guidelines (Appendix 2 of the
Nano Letters author guide), centralized in `common/acs_style.py`:

- Single-panel figures are single-column width (exactly 3.33 in / 240 pt, the
  ACS maximum); any figure with more than one panel is double-column width
  (exactly 7.0 in / 504 pt, likewise the maximum). Depths are well inside the
  9.167 in cap.
- One typeface and one size everywhere: 8 pt Arial, for text and for mathtext,
  against an ACS floor of 4.5 pt.
- Nothing is drawn thinner than the 0.5 pt ACS floor: 0.5 pt spines, ticks,
  grids and error bars; 0.6 pt marker edges; 1.0 pt data and fit curves.
  Markers are 3.5 pt, 4.5 pt where marker shape encodes a variable, and 2.4 pt
  (with a 0.5 pt edge) on the dense profiles that carry several hundred points
  per curve.
- Series are separated by marker shape (or a colourbar) as well as colour, so
  no information is carried by colour alone.
- Fonts are embedded as TrueType (`pdf.fonttype = 42`), never Type 3.
- `save_fig()` writes the PDF at exactly the figure size, with no
  `bbox_inches='tight'` trimming, so the media box stays exactly one or two
  columns wide and the manuscript can embed it unscaled — 8 pt drawn is 8 pt
  printed.
- Multi-panel figures are labeled (a), (b), (c)... via `label_panels()`.

Regenerate a figure by running its script directly from its own directory,
e.g. `cd SAXPCS_8id && python3 saxpcs.py` (headless environments can set
`MPLBACKEND=Agg` to skip the interactive `plt.show()`).
