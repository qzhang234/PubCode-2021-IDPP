# Dynamic Arrest Follows Assembly on a Temperature-Scaled Common Clock in an Elastin-like Polypeptide

Everything needed to reproduce the PNAS submission: the reduced experimental
data, the analysis code that turns it into every data-derived figure, and the
LaTeX sources for the manuscript, the SI Appendix and the cover letter.

The measurement is small-angle X-ray photon correlation spectroscopy (SA-XPCS)
at APS beamline **8-ID-I**, October 2021 (beamtime 2021-3): four loadings of
the elastin-like polypeptide (VPAVG)30 held at 28.00, 28.25, 28.50 and
28.75 °C. Correlating the detector frames gives the fast and arrested
fractions at five length scales; averaging the same frames gives I(Q).

```
make env        # create the conda environment, once
make figures    # regenerate every data-derived figure from the reduced data
make papers     # compile main.pdf, si.pdf, cover_letter.pdf into manuscript/build/
make all        # figures + papers
make check      # report undefined citations / references / missing figures
make robustness # recompute the robustness numbers of SI sections 4.3 and 6.2 (no figure)
make clean      # remove LaTeX aux files, keep the PDFs
```

## Which code makes which figure

The figure files keep the names of an earlier figure order. The manuscript
cites them by label, so the printed number is the one in the first column.
`make figures` runs every script below from the committed data. `make papers`
then finds each PDF where its script wrote it, because the Makefile puts those
directories on `TEXINPUTS`.

| Printed figure | File | Made by | Reads |
|---|---|---|---|
| Fig. 1 | `manuscript/figures/Setup.pdf` | not generated: the approved illustrator schematic | -- |
| Fig. 2 | `analysis/SAXPCS_8id/Figure3_Isothermal_SAXPCS.pdf` | `analysis/SAXPCS_8id/temperature_series_2021-3.py` | `data_2021-3/` |
| Fig. 3 | `analysis/SAXPCS_8id/Figure4_Arrest_vs_Temperature.pdf` | `analysis/SAXPCS_8id/arrest_vs_temperature.py` | `data_2021-3/` |
| Fig. S1 | (a) typeset in `si.tex`; (b) `manuscript/figures/FigureS1b_SDS_PAGE.png` | not generated: primer and template sequences, and a gel photograph | -- |
| Fig. S2 | `analysis/SAXPCS_8id/FigureS2_Calibration.pdf` | `analysis/SAXPCS_8id/calibration_figure.py` | constants in `abs_xsec.py`; no data file |
| Fig. S3 | `analysis/Rad_Dam_Check/FigureS6_Flux_Control.pdf` | `analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py` | `Rad_Dam_Check/cluster_results/`, `Pind_calibration.csv` |
| Figs. S4-S7 | `analysis/SAXPCS_8id/FigureS7_Temperature_2800.pdf` to `FigureS10_Temperature_2875.pdf`, one per temperature | `temperature_series_2021-3.py` | `data_2021-3/` |
| Fig. S8 | `analysis/SAXPCS_8id/FigureS12_Fast_Mode.pdf` | `temperature_series_2021-3.py` | `data_2021-3/` |
| Fig. S9 | `analysis/SAXPCS_8id/FigureS11_Collapse_AllQ.pdf` | `arrest_vs_temperature.py` | `data_2021-3/` |

Fig. 2 is the top row of Fig. S6: the same 28.50 °C fits, drawn twice by one
script. Every correlation-function fit behind Figs. 2, 3 and S4-S9 is Equation
2 of the main text, fitted by one function, `fit_seeded()` in
`analysis/SAXPCS_8id/xpcs_fit.py`, in both figure scripts.

Some numbers in the text come from scripts that print them rather than draw:

| Script | Prints | Run by |
|---|---|---|
| `analysis/SAXPCS_8id/temperature_series_2021-3.py` | the numbers quoted with Fig. 2 and in SI section 6: fast fractions, D0, low-Q intensities of every group, background checks | `make figures` |
| `analysis/SAXPCS_8id/arrest_vs_temperature.py` | T_s and the master curve, the collapse at each Q, the SAXS on the same clock, the half-arrest intensities, and the checks on the half-arrest threshold and the Q band; then E_a = k_B T^2/T_s, the Eyring free energy and entropy of activation from the half-arrest times, D0, the fast relaxation time 1/(D0 Q^2) and the waiting time t_w compared at the same rescaled time across temperatures, and the bond number N from published calorimetric enthalpies, scaled to 30 repeats with the chain-length dependence read from Girotti et al. 2004, Fig. 4 (SI section 6.6) | `make figures` |
| `analysis/SAXPCS_8id/general_model_check.py` | how far T_s and f move under five versions of the g2 fit (SI sections 4.3 and 6.2); about six minutes | `make robustness` |
| `analysis/SAXPCS_8id/contrast_2021-3.py` | the speckle contrast, beta = 0.141, from the static standard | by hand |
| `analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py` | the fluxes and ratios quoted with Fig. S3 | `make figures` |
| `analysis/Rad_Dam_Check/Flux_Cal.py` | the raw flux for the recorded PIN-diode reading | by hand |

## Code structure

```
.
├── Makefile              the whole build: figures, the three PDFs, make check, make robustness
├── README.md             this file
├── analysis/
│   ├── environment.yml   the conda environment env-2021-LLDP, matplotlib pinned
│   ├── common/
│   │   ├── acs_style.py  house figure style: column widths, 8 pt Arial, line widths, save_fig()
│   │   └── utils.py      the first outlier cut and the group average, used by the reduction
│   ├── SAXPCS_8id/       the 2021-3 temperature series
│   │   ├── data_2021-3/                  committed reduced data (127 MB): 85 averaged NeXus files
│   │   │                                 and trace_2021-3.csv, one row per acquisition
│   │   ├── temperature_series_2021-3.py  figure script: Fig. 2, Figs. S4-S8
│   │   ├── arrest_vs_temperature.py      figure script: Fig. 3, Fig. S9
│   │   ├── calibration_figure.py         figure script: Fig. S2
│   │   ├── general_model_check.py        robustness numbers, no figure
│   │   ├── contrast_2021-3.py            speckle contrast, no figure
│   │   ├── xpcs_fit.py       shared: Equation 2, the seeded fit, beta from the static standard
│   │   ├── series_2021_3.py  shared: the t_w clock and the absolute, background-subtracted I(Q)
│   │   ├── nexus_read.py     shared: readers for the averaged NeXus files
│   │   ├── abs_xsec.py       shared: the ion-chamber calibration to absolute intensity (mm^-1)
│   │   ├── reduce_babnigg202110.py, nexus_convert_babnigg202110.py,
│   │   │   restore_dx_metadata_2021-3.py, make_qmap_babnigg202110.py,
│   │   │   sweep_stuck_pixels_2021-3.py, run_ and stop_babnigg202110_correlations.sh
│   │   │                     beamline-side provenance, never run by make (see below)
│   │   ├── *_manifest.csv, nexus_manual.txt   records of the NeXus conversion
│   │   └── Figure*.pdf       the generated figures
│   └── Rad_Dam_Check/        the flux-dependence control; see its README.md
│       ├── g2_SAXPCS_Rad_Cali.py   figure script: Fig. S3
│       ├── Flux_Cal.py             PIN-diode current to photon flux
│       ├── Pind_calibration.csv, PINdiode_Calibration_DP00429_2024-02-01.pdf
│       └── cluster_results/        the four committed result files, E0171-E0174
├── manuscript/           LaTeX sources and compiled PDFs; see manuscript/README.md
│   ├── main.tex, si.tex, cover_letter.tex, cover_letter_header.tex
│   ├── reference.bib
│   ├── figures/          artwork that is not generated: Setup.pdf (Fig. 1),
│   │                     FigureS1b_SDS_PAGE.png (Fig. S1b), anl.png (letterhead logo);
│   │                     the other three images are not used by the current documents
│   ├── build/            main.pdf, si.pdf and cover_letter.pdf, tracked; aux files, ignored
│   └── texmf/            vendored achemso and captdef, so the build needs no root install
└── tools/
    ├── check_refs.py     checks every cited reference against CrossRef
    └── pdf_comments.py   matches exported PDF review comments to the text they mark
```

`data_2021-3/` holds the group averages of the four loadings and their 6 °C
references, the empty-cell and air measurements, and the static standard used
for the speckle contrast.

Kept locally and gitignored: `PNAS_author_guide.pdf`, `ACS_author_guide.pdf`
(Appendix 2 of the ACS guide is the source of the numbers in
`analysis/common/acs_style.py`), and `References/`, the PDFs of cited papers
used for the citation audit. `.loglogin` and `analysis/.loglogin` are beamline
login logs carried over from the acquisition account; nothing here uses them.

### Data flow

```
/gdata/s8id-dmdtn/2021-3/babnigg202110         raw frames and 2021 metadata (beamline only)
    --> nexus_convert_babnigg202110.py, restore_dx_metadata_2021-3.py
/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus   NeXus metadata
    --> make_qmap_babnigg202110.py, sweep_stuck_pixels_2021-3.py     the q-map and its mask
    --> boost_corr_bin (run_babnigg202110_correlations.sh)          one result per acquisition
    --> reduce_babnigg202110.py      groups of 100 acquisitions, two outlier cuts
analysis/SAXPCS_8id/data_2021-3/               committed; everything below runs from here
    --> temperature_series_2021-3.py  --> Fig. 2, Figs. S4-S8
    --> arrest_vs_temperature.py      --> Fig. 3, Fig. S9
    --> general_model_check.py, contrast_2021-3.py  --> numbers on stdout
abs_xsec.py constants                 --> calibration_figure.py   --> Fig. S2
analysis/Rad_Dam_Check/cluster_results/ --> g2_SAXPCS_Rad_Cali.py --> Fig. S3
manuscript/*.tex and the figures above --> make papers --> manuscript/build/*.pdf
```

The figure scripts share their building blocks instead of copying them:
`xpcs_fit.py` holds the model and the fit, `series_2021_3.py` the clock and
the SAXS, `nexus_read.py` the file readers, and `abs_xsec.py` the absolute
scale. Every figure script imports `common/acs_style.py`, so all the figures
share one style. Each script is a plain top-to-bottom program; run it from any
directory, and it reads and writes next to itself.

## How the figures are computed

The analysis choices, and where each lives in the code. The SI Appendix gives
the evidence for them; the numbers quoted here are ones the scripts print.

**The waiting-time clock** (`series_2021_3.py`). t_w = 0 is the first
acquisition whose stage-RTD reading, as recorded in `trace_2021-3.csv`, is
within 0.01 °C of its target. Acquisitions follow one another every 6 s, each
at a fresh sample position, and are averaged in groups of 100, with one group
boundary placed at the first acquisition on target. A group is placed at the
midpoint of the acquisitions it averages. Placing it at its first acquisition
instead would put the first group at about 6 s rather than about 300 s, more
than a decade early on the logarithmic axes of Fig. 3, and would bias T_s low. Groups
that start before t_w = 0 belong to the ramp and are not fitted.

**Absolute intensity** (`abs_xsec.py`; SI section 3). I_abs(Q) = N_frames /
(t_exp F T L dOmega) x I(Q), following Sheyfer et al., Phys. Rev. Lett. 125,
125504 (2020), Eq. S1, with I(Q) in photons per pixel per frame, t_exp = 2 s,
L = 2 mm and dOmega = (76 um / 7.8 m)^2. Each averaged file carries the mean
upstream and downstream ion-chamber readings of its acquisitions, so every
curve is scaled by its own flux F and transmission T. F follows from the
upstream chamber through a PIN-diode calibration, a straight line fitted to
seven of eight attenuation points (the unattenuated point lies off the line
the others define); T is the ratio of the two chambers divided by the
transmission of the air path between them, 0.8679 +/- 0.0009. The calibration
was recorded in a separate beamtime on the same ion chambers; its eight points
are constants of `abs_xsec.py` and are drawn in Fig. S2. Lengths are in
millimetres, so every absolute axis is in mm^-1, ten times smaller than the
cm^-1 the literature more often quotes. The implementation reproduces the
original beamline notebooks (not in this repository) to 2% when the thickness
is set to their 1 mm.

**Background subtraction** (`series_2021_3.py`, `temperature_series_2021-3.py`;
SI section 6.2). The empty cell E0411 is subtracted on the absolute scale with
the coefficient `BG_COEF` = 0.80. Six candidate background measurements exist,
and their beamline labels do not reliably say what was in the cell, so the
choice rests on the monitors and on I(Q): E0411 is the only one with material
in the beam that the chambers see, and none holds a liquid column of the
sample thickness. The coefficient trades the flatness of the 6 °C references
against a subtraction dip near 0.004-0.005 A^-1. To compare coefficients, each
reference is reduced to its median in nine logarithmic Q bins over
0.0036-0.0305 A^-1 and divided by its plateau over 0.012-0.0305 A^-1, so every
decade of Q counts equally; least squares on I itself would be dominated by
the lowest Q, where the empty cell is several times the 6 °C signal. The
flattest coefficient is 0.847-0.882 for the four references; at 0.80 they rise
1.26- to 1.45-fold toward low Q and turn back up by at most 4%. No buffer was
measured in this beamtime, so the subtraction cannot remove windows, solvent
and parasitic scattering together: at 0.004 A^-1 the subtracted background is
2.4-3.0 times the 6 °C signal. The fitted groups average 54-322 times the
subtracted background over 0.004-0.008 A^-1, so this limit bears on the 6 °C
references alone. I(Q) is drawn from the lowest fitted Q bin outward; below it
the points sit on the edge of the beamstop shadow. Loading D8 was measured at
6 °C before and after its isothermal experiment, and the earlier measurement
(H0402, the first of the beamtime) still holds undissolved material, so D8's
reference is the later one, H0403.

**Speckle contrast** (`xpcs_fit.measure_beta()`, `contrast_2021-3.py`; SI
section 6.2). beta is a property of the beam and the detector, so it is
measured once and held fixed in every fit. The correlation function of a
static sample is flat at 1 + beta. On the static standard B0005 (NIST glassy
carbon) the flat level is taken in every Q bin, and beta = 0.141 is the
inverse-variance mean over the nine bins between 0.0085 and 0.0165 A^-1. The
five fitted bins read lower, 0.1226-0.1390: the standard is scratched, and its
static parasitic scattering depresses the apparent contrast most at low Q. The
band stops before a dip at 0.0175-0.0185 A^-1 and before the bin of a stuck
pixel at 0.0214 A^-1. One scalar serves every bin and every group, because
inside the plateau the bin-to-bin scatter is noise that a per-bin beta would
feed into f.

**The correlation-function fit** (`xpcs_fit.py`; SI sections 4.1-4.3 and 6.3).
Equation 2 is fitted to the five lowest dynamic Q bins, 0.00359-0.00756 A^-1
(175-83 nm); above them the intensity has fallen by more than an order of
magnitude and g2 is dominated by counting noise. Within a group, f and
tau_slow are free in each bin, while D0, p_fast and p_slow are shared by the
five bins; groups are fitted independently. The fast time is tied to
tau_fast = 1/(D0 Q^2), a form imposed rather than fitted. At that delay the
fast term has fallen to 1/e of its amplitude whatever p_fast is, so D0 is the
rate of the fast decay and compares groups whose stretching differs (a mean
relaxation time would carry a further factor Gamma(1/p)/p, 1 at p = 1 but 120
at p = 0.2, set by the tail of the decay). The fit
minimises the error-weighted residuals with `scipy.optimize.least_squares`,
using the stored g2 uncertainties, within the bounds 0 <= f <= 1,
1 <= tau_slow <= 1e4 s, 0.2 <= p <= 3 and 1e2 <= D0 <= 1e10 A^2/s. Parameter
uncertainties are the square roots of the diagonal of the pseudo-inverse of
J^T J, so a parameter that settles on a bound (f = 1 in some mobile groups)
leaves the others defined. With p_fast free the fit has local minima in which
the fast term imitates the slow one, so `fit_seeded()` first fits with p_fast
held at each of 0.20, 0.25, ..., 1.30, then releases it from the best of those
and keeps the released fit only if its chi^2 is no worse. p_fast is free in
each group, because it falls from about 1 to about 0.2 as a sample arrests.

**Which groups are fitted, and what is drawn** (`temperature_series_2021-3.py`).
A group is fitted when the median g2 uncertainty at the lowest Q between 1e-4
and 1e-2 s is at most 0.016, a noise limit carried over from an earlier
analysis, and kept when its reduced chi^2 is below 10; 34 groups pass. Where f
at the lowest Q is below 0.25 the fast decay is too weak to fix a rate, and D0
is drawn as an open grey star and left out of the quoted rate ranges. Arrest
counts as under way once f at the lowest Q falls below 0.85. When an
isothermal experiment has more than seven fitted groups, panels (a), (b) and
(d)-(g) show every second one and panel (c) all of them; the group at t_w = 0
is drawn black in panel (a) when it is too noisy to fit.

**The temperature collapse and the SAXS on the same clock**
(`arrest_vs_temperature.py`; SI section 6.5). f at the lowest Q of all 34
groups is fitted with one master curve, f = A / [1 + (X/X_c)^n] with
X = t_w exp[(T - T0)/T_s], each point weighted by its fit uncertainty floored
at 0.01. T0 = 28.00 °C is a convention: exp[(T - T0)/T_s] factorises, so T0
slides the curve along X and leaves T_s unchanged. The uncertainty of T_s is
the standard deviation over 200 bootstrap resamples of whole groups (random
seed 0). Each higher Q bin is collapsed on its own the same way. The low-Q
intensity is the mean absolute I(Q) of each group over 0.004-0.008 A^-1, the
length scales of the fitted bins, where I(Q) grows most; every group at the
target is used, fitted or not. The T_s that superposes the intensity on its
own minimises, on a grid from 0.10 to 1.50 K, the scatter of log I about one
cubic in log X. A sample is half arrested when f at the lowest Q falls to 0.5,
interpolated linearly between the two groups that bracket it. The script also
prints the same reading for thresholds of 0.3 to 0.7 and for four other Q
bands.

**Robustness** (`general_model_check.py`, `make robustness`). The collapse is
repeated with f from five versions of the fit, listed in the script's
docstring: Equation 2 as used; the general two-mode form with tau_fast free in
every bin, seeded and unseeded; p_fast fixed at 0.5; and the lower contrast of
the fitted bins.

## Reproducibility

The reduced data every figure is drawn from are committed:
`analysis/SAXPCS_8id/data_2021-3/` (127 MB: one averaged NeXus file per group
of acquisitions, plus the per-acquisition trace of time, RTD temperature,
monitors and low-Q intensity) and `analysis/Rad_Dam_Check/cluster_results/`.
**`make figures` and `make papers` therefore run from a clone of this
repository alone** -- no beamline account, no network, no raw data.

### Code that reads outside this repository

These scripts document how the committed data were produced. None of them is
on the path from the committed data to a figure, and `make figures` never runs
them. Every beamline data path is on `/gdata`, which APS backs up and which
mounts only on the beamline GPU hosts (amber, amazonite, adamite). They need
packages beyond `environment.yml`: `nexus_xpcs_aps`, the pySimpleMask source,
and `boost_corr`.

| Script | Out-of-repo location |
|---|---|
| `analysis/SAXPCS_8id/nexus_convert_babnigg202110.py` | `/gdata/s8id-dmdtn/2021-3/babnigg202110` (the Data Exchange archive) and `..._nexus` (the NeXus tree it writes) |
| `analysis/SAXPCS_8id/restore_dx_metadata_2021-3.py` | the same pair; checks that every file carries its monitors |
| `analysis/SAXPCS_8id/make_qmap_babnigg202110.py` | the same pair, the detector blemish map and the pySimpleMask source under `/home/beams*/8IDIUSER` |
| `analysis/SAXPCS_8id/sweep_stuck_pixels_2021-3.py` | the NeXus tree |
| `analysis/SAXPCS_8id/run_babnigg202110_correlations.sh`, `stop_babnigg202110_correlations.sh` | the NeXus tree; drive `boost_corr_bin` on a GPU host |
| `analysis/SAXPCS_8id/reduce_babnigg202110.py` | `/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results`, the per-acquisition results `data_2021-3/` is averaged from. It also imports `average_ranges.py`, which is not part of this repository and which uses `analysis/common/utils.py`. |

All paths are hard-coded at the top of each script. Every other script anchors
its input and output on its own directory and runs from any working directory.

### How the committed data were produced

**NeXus conversion** (`nexus_convert_babnigg202110.py`, manifests alongside).
8-ID-I moved to NeXus around 2024 and `boost_corr_bin` reads only NeXus
metadata, so the 2021 Data Exchange records were converted. The frames are
symlinked, not copied. The ion chambers, ring current, times, geometry and
sample temperatures are read from the archive, and a field an acquisition did
not record is written as NaN, never as a placeholder. A few acquisitions
without a complete Data Exchange record are listed in the manifests and were
not converted. Data Exchange indexes the detector with `beam_center_x` along
its 512-pixel axis and NeXus with x along the 1024-pixel axis, so the
archive's centre (x=260.6, y=512.2) becomes (x=512.2, y=260.6); refined in
pySimpleMask against the measured scattering, it is (x=510.5, y=259.5).
`beam_center_position_x/y` are set equal to `position_x/y`, so the 8-ID-I
reader adds no offset to the centre. `restore_dx_metadata_2021-3.py` checks
that every metadata and result file carries the archive's monitors.

**The q-map** (`make_qmap_babnigg202110.py`). Geometry from the archive:
10.945 keV, 7.8 m, 76 um pixels, 512 x 1024. Mask: the union of the
hand-drawn 2021-3 beamstop, the detector blemish map, the pixels the 2021
Matlab mask excluded, and two stuck pixels, then q >= 0.0032 A^-1. Partition:
270 x 8 static and 27 x 1 dynamic bins, linear. The series is correlated
against `babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm2.hdf`.
CrossRm was built in the pySimpleMask 0.4.0 GUI, which adds the beamstop
cross arms and refines the beam centre; CrossRm2 is CrossRm with the second
stuck pixel masked, made by copying CrossRm and clearing that one pixel, so
every q value and bin edge is unchanged. The first stuck pixel, (30, 340),
reads about 0.75 photons per frame in every measurement. The second, (63, 166),
reads 1.5 photons per frame in every measurement from 17 Oct on, lies at
Q = 0.0214 A^-1, outside the five fitted bins, and was found by
`sweep_stuck_pixels_2021-3.py`, which compares every pixel with the median of
its q ring in every scan prefix of the beamtime.

**Correlation** (`run_babnigg202110_correlations.sh`). `boost_corr_bin`,
multi-tau, one GPU per process, with `--crop-ratio-threshold 0.5` (the default
of 1.0 crops whenever any pixel is masked, which returns all-NaN g2 for these
data) and `-w`, which overwrites a result in place.

**Reduction** (`reduce_babnigg202110.py`). Writes `trace_2021-3.csv`, one row
per acquisition with its start time, zone-3 RTD reading and setpoint, both ion
chambers and the mean of the stored I(Q) entries over 0.004-0.008 A^-1; the
groups of 100 acquisitions; each 6 °C reference and the static standard,
averaged whole; and, with the argument `empties`, the six candidate
backgrounds. Run both calls after a reprocessing. Every group passes two
outlier cuts (SI section 5.3). The first, `outlier_removal()` in
`common/utils.py`, compares the shape of each acquisition's g2, g2
uncertainty and log I(Q) with the rest of the group by cosine similarity and
drops the least similar 5% on each. The second, in `average_ranges.py`, drops
an acquisition whose mean I(Q) over 0.004-0.008 A^-1 exceeds the group median
by more than three robust standard deviations: a transient scatterer that the
shape score misses.

### Where the 8-ID-I data live

| Path | What it holds |
|---|---|
| `/gdata/s8id-dmdtn/2021-3/babnigg202110` | Raw detector frames and the 2021 Data Exchange metadata. The archival copy, never written to. |
| `/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus` | NeXus metadata rewritten for the current `boost_corr_bin`, the q-map, and `reprocess_results/`, the per-acquisition correlation results. The `.bin` files are symlinks into the archive. |
| `/gdata/s8id-dmdtn/2025-1/qzhang202503/flux_control_results` | The four flux-control result files, E0171–E0174 (2025). Their raw frames are in the DM experiment qzhang202503, archived on 23 Aug 2025; see `analysis/Rad_Dam_Check/README.md`. |

## Source data integrity

No data file is modified by producing a figure. The figure and number scripts
and the modules they import open every data file read-only. The scripts that
write (`make_qmap_babnigg202110.py`, `nexus_convert_babnigg202110.py`,
`restore_dx_metadata_2021-3.py --apply` and `reduce_babnigg202110.py`) write
to beamline storage or rebuild `data_2021-3/`, and `make` never runs them.

Nothing is rescaled between the calibration and the plot. `abs_xsec_coef()` is
built from millimetres, so every absolute intensity axis in this work is
**mm⁻¹**, plotted as the coefficient comes out. (The literature more often
quotes cm⁻¹, which is ten times larger.)

## Submission target

PNAS Research Report, Direct Submission (`PNAS_author_guide.pdf`): preferred
length six pages, about 4,000 words, 50 references and four medium figures;
abstract at most 250 words; significance statement at most 120 words; at least
three keywords; SI Appendix as a single PDF.
