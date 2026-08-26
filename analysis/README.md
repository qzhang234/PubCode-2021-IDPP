# llps-saxpcs

Data reduction and figure-generation code for the manuscript "Reversible
β-Sheet-like Ordering and XPCS-Resolved Dynamic Arrest in Phase-Separating
Elastin-like Polypeptide Networks" (main text + Supplemental Information).

The experiment combines small/wide-angle X-ray scattering (SAXS/WAXS) at
12-ID-B with X-ray photon correlation spectroscopy (XPCS) at 8-ID-I to follow
an elastin-like polypeptide (ELP) sample through a liquid-liquid phase
separation (LLPS): SAXS/WAXS resolves the static structure (Guinier radius,
the emergent β-sheet-like correlation peaks), while XPCS resolves the
dynamics (g2 relaxation) as the sample ages after a temperature jump. This
repo goes from raw beamline files to the PDF figures used in the paper.

## Repository layout

- **`common/`** — code shared by more than one plotting script.
  - `utils.py`: parallel HDF reading, cross-correlation-based outlier removal,
    and averaging helpers for XPCS cluster-result files.
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
raw .avg files (12-ID-B, /home/8-id-i/2021-1/12-id-b/.../Processed/)
    --> Read_12ID_SAWAXS.py   (merge SAXS+WAXS, scale WAXS onto SAXS, propagate errors)
    --> reduced_data/Merged_*.csv
    --> Guinier_Plot.py       --> FigureS7_Guinier.pdf
    --> Plot_12ID.py          --> Figure2_SAXS_WAXS.pdf
```

`Guinier_Plot.py` fits `ln[I(Q)]` vs. `Q^2` in the low-Q Guinier region of the
10 °C reference to get R_g. `Plot_12ID.py` overlays the merged reference and
measurement curves, annotates the high-Q power law, and fits the two
β-sheet-like correlation peaks with Gaussians in an inset.

### `SAXPCS_8id/`

```
raw cluster-result .hdf files (8-ID-I, per-file g2, g2_err, saxs_1d, IC readout)
    --> saxpcs.py            --> Figure3_Isothermal_SAXPCS.pdf,
                                  FigureS9_Fit_Parameters.pdf,
                                  FigureS3_Calibration.pdf
    --> saxs_evolution.py    --> FigureS8_SAXS_Evolution.pdf
    --> g2_grid_SI.py        --> FigureS10_g2_Grid.pdf
```

`saxpcs.py` reads the averaged B0147 (isothermal 30 °C) and D0138 (buffer)
HDF files directly with h5py, fits a double-exponential (KWW) model to g2(τ)
at the primary q bin, then fits each elapsed time's τ_fast(Q) and τ_slow(Q)
to a power law τ = A·Q^γ (weighted log-log regression with error propagation)
to extract the scaling exponents γ_fast and γ_slow. It also converts SAXS
intensities to an absolute differential cross section via `abs_xsec.py`
(per-file ion-chamber-based coefficients — see the calibration derivation in
that file), and calibrates the upstream ion-chamber reading against the
incident photon flux. `saxs_evolution.py`
reuses the same absolute-cross-section calibration to show the full SAXS
evolution across all B0147 files. `g2_grid_SI.py` repeats the g2 fit for the
four secondary q bins not shown in the main figure. `average_ranges.py` is a
helper for building custom-range averages instead of the default per-file
grouping.

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
artifact. `Flux_Cal.py` holds the underlying flux-conversion formula.

## Figure-to-manuscript mapping

Each output PDF is named exactly as `main.tex` / `si.tex` reference it, so the
LaTeX build picks the files up straight from this directory — only the path
differs, and the top-level `Makefile` resolves it via `TEXINPUTS`.

| Output PDF | Script | Figure | Panels |
|---|---|---|---|
| `SAXS_12id/Figure2_SAXS_WAXS.pdf` | `Plot_12ID.py` | Main Fig. 2 | 1 (+ inset) |
| `SAXS_12id/FigureS7_Guinier.pdf` | `Guinier_Plot.py` | SI Fig. S5 | 1 |
| `SAXPCS_8id/Figure3_Isothermal_SAXPCS.pdf` | `saxpcs.py` | Main Fig. 3 | 3 |
| `SAXPCS_8id/FigureS9_Fit_Parameters.pdf` | `saxpcs.py` | SI Fig. S7 | 4 (2×2) |
| `SAXPCS_8id/FigureS3_Calibration.pdf` | `saxpcs.py` | SI Fig. S3 | 2 |
| `SAXPCS_8id/FigureS8_SAXS_Evolution.pdf` | `saxs_evolution.py` | SI Fig. S6 | 1 |
| `SAXPCS_8id/FigureS10_g2_Grid.pdf` | `g2_grid_SI.py` | SI Fig. S8 | 4 (2×2) |
| `Rad_Dam_Check/FigureS5_Flux_Control.pdf` | `g2_SAXPCS_Rad_Cali.py` | SI Fig. S4 | 2 |

Main Fig. 1 (`manuscript/figures/Setup.pdf`) is an illustrator schematic and is
not produced by any script here. SI Fig. S1 and the Video S1 still frame are
photographs/renderings, also in `manuscript/figures/`. SI Fig. S2 is a
LaTeX-typeset sequence box with no image file.

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
  Markers are 3.5 pt, or 4.5 pt where marker shape encodes a variable.
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
