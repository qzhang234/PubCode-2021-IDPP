# analysis

Data reduction and figure generation for "Dynamic Arrest Follows Assembly on a
Temperature-Scaled Common Clock in an Elastin-like Polypeptide" (main text and
SI Appendix). Everything here comes from one SA-XPCS beamtime at APS 8-ID-I,
October 2021 (2021-3), apart from the flux-dependence control in
`Rad_Dam_Check/`.

## Layout

- **`common/`** -- `acs_style.py`, the figure-formatting module every plotting
  script imports (column widths, 8 pt Arial, line and marker sizes, panel
  letters, `save_fig()`); `utils.py`, the outlier cuts and averaging used by
  the reduction, plus beamline-side helpers that cannot run from a clone.
- **`SAXPCS_8id/`** -- the 2021-3 temperature series, from the raw archive to
  the figures.
- **`Rad_Dam_Check/`** -- the flux-dependence control. See its README.

## Data flow

```
2021-3 Data Exchange archive (/gdata/s8id-dmdtn/2021-3/babnigg202110; BEAMLINE-ONLY)
    --> make_qmap_babnigg202110.py        geometry and mask provenance;
                                          --crossrm2 writes the q-map in use
    --> sweep_stuck_pixels_2021-3.py      finds the stuck pixels that mask needs
    --> nexus_convert_babnigg202110.py    NeXus tree, monitors carried across
    --> restore_dx_metadata_2021-3.py     check that every file carries them
    --> boost_corr_bin                    (run_babnigg202110_correlations.sh)
    --> reduce_babnigg202110.py           --> data_2021-3/   (groups of 100 acquisitions,
                                                              two outlier cuts, and the trace)
data_2021-3/
    --> temperature_series_2021-3.py  --> Figure3_Isothermal_SAXPCS.pdf (main Fig. 2),
                                          FigureS7-S10_Temperature_*.pdf (SI Figs. S4-S7),
                                          FigureS12_Fast_Mode.pdf (SI Fig. S8)
    --> arrest_vs_temperature.py      --> Figure4_Arrest_vs_Temperature.pdf (main Fig. 3),
                                          FigureS11_Collapse_AllQ.pdf (SI Fig. S9)
    --> general_model_check.py        --> robustness numbers of SI sections 4.3 and 6.2
                                          on stdout (no figure)
    --> contrast_2021-3.py            --> the speckle contrast on stdout (no figure)
abs_xsec.py constants
    --> calibration_figure.py         --> FigureS2_Calibration.pdf (SI Fig. S2)
Rad_Dam_Check/cluster_results/
    --> Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py --> FigureS6_Flux_Control.pdf (SI Fig. S3)
```

The file names carry the numbers of an earlier figure order; the manuscript
resolves them by label, so the names are not the printed numbers. The printed
number of each figure is given in parentheses above.

## The shared modules

- `xpcs_fit.py` -- the two-mode g2 model of Equation 2 and its global fit, and
  `measure_beta()`, which reads the speckle contrast off the plateau of the
  static standard (the standard's scratches depress it in the fitted bins).
- `abs_xsec.py` -- the ion-chamber calibration behind every absolute intensity.
- `series_2021_3.py` -- the clock and the SAXS shared by the two figure
  scripts: t_w = 0 at the first acquisition with the RTD at its target, each
  group placed at the midpoint of the acquisitions it averages, and the
  absolute I(Q) with the empty cell E0411 subtracted at the coefficient its
  comments justify.
- `nexus_read.py` -- reading the averaged NeXus files.

## What the figure scripts compute

`temperature_series_2021-3.py` fits Equation 2 to the five lowest Q bins of
every group, with beta fixed, p_fast free and started from the minimum of a
p_fast profile. It draws the four per-temperature figures, SI Figs. S4-S7
(the 28.50 C one is also main Fig. 2), and the fast mode of every group
against the arrested fraction, SI Fig. S8. It prints the numbers SI section 6
quotes. Its docstring lists
every analysis choice and why, including the background measurement and why
loading D8 has two 6 C references.

`arrest_vs_temperature.py` refits the same groups and collapses the fast
fraction at the lowest Q of all four holds onto one master curve,
f = A / [1 + (X/X_c)^n] with X = t_w exp[(T - T0)/T_s]. It bootstraps T_s over
whole groups, collapses each higher Q bin on its own for SI Fig. S9, places the
low-Q SAXS of the same groups on the same clock, and prints the numbers main
Figure 3 and the scaling paragraphs quote, including the checks that the
half-arrest threshold and the choice of Q band do not change the reading of
Figure 3c.

`general_model_check.py` repeats that collapse with f from five versions of
the g2 fit (Equation 2; the general two-mode form seeded and unseeded;
p_fast frozen at 0.5; the lower contrast of the fitted bins). Its numbers are
quoted in SI sections 4.3 and 6.2. It takes about a quarter of an hour and is
run by `make robustness`, not by `make figures`.

## Environment

`environment.yml` pins the scientific packages the scripts need:

```
conda env create -f environment.yml     # creates env-2021-LLDP
```

matplotlib is pinned to 3.9.4 so layout does not move between releases.
`pyxpcsviewer` is not required: it is imported lazily inside the beamline-side
helpers of `common/utils.py`, which no figure depends on.

## Figure formatting

`common/acs_style.py` fixes the house style for every figure: single-column
figures exactly 3.33 in wide and double-column exactly 7.0 in, 8 pt Arial for
every character, nothing thinner than 0.5 pt, fonts embedded as TrueType, and
series separated by marker shape as well as colour. `save_fig()` writes at the
exact figure size, so the manuscript embeds each figure unscaled and 8 pt drawn
is 8 pt printed. These sizes also satisfy PNAS (text 6-12 pt; figures up to
18 cm wide).
