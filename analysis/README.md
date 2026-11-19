# analysis

Data reduction and figure generation for the paper. The repository-root
`README.md` is the map: which script makes which figure, what each one reads
and writes, the order to run them, and how every figure is computed.

- **`common/`** -- `acs_style.py`, the figure style every plotting script
  imports; `utils.py`, the first outlier cut and the group average used by the
  beamline-side reduction.
- **`SAXPCS_8id/`** -- the 2021-3 temperature series: the committed reduced
  data (`data_2021-3/`), the figure scripts and the modules they share, and the
  beamline-side scripts that produced the data.
- **`Rad_Dam_Check/`** -- the flux-dependence control, SI Fig. S3; see its
  README.
- **`environment.yml`** -- the conda environment of every figure and number
  script (`make env` from the repository root). matplotlib is pinned to 3.9.4
  so the layout of the figures does not move between releases.

The figures follow the house style set in `common/acs_style.py`: single-column
figures exactly 3.33 in wide and double-column exactly 7.0 in, 8 pt Arial for
every character, nothing thinner than 0.5 pt, and fonts embedded as TrueType.
`save_fig()` writes at the exact figure size, so the manuscript embeds each
figure unscaled and 8 pt drawn is 8 pt printed.
