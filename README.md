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
2 of the main text. It is fitted by `fit_g2_diffusive()` in
`analysis/SAXPCS_8id/xpcs_fit.py`, called through the same seeding wrapper,
`fit_block()`, in both figure scripts.

Some numbers in the text come from scripts that print them rather than draw:

| Script | Prints | Run by |
|---|---|---|
| `analysis/SAXPCS_8id/temperature_series_2021-3.py` | the numbers quoted with Fig. 2 and in SI section 6: fast fractions, D0, low-Q intensities of every group, background checks | `make figures` |
| `analysis/SAXPCS_8id/arrest_vs_temperature.py` | T_s and the master curve, the collapse at each Q, the SAXS on the same clock, the half-arrest intensities, and the checks on the half-arrest threshold and the Q band | `make figures` |
| `analysis/SAXPCS_8id/general_model_check.py` | how far T_s moves under five versions of the g2 fit (SI sections 4.3 and 6.2); about a quarter of an hour | `make robustness` |
| `analysis/SAXPCS_8id/contrast_2021-3.py` | the speckle contrast, beta = 0.141, from the static standard | by hand |
| `analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py` | the fluxes and ratios quoted with Fig. S3 | `make figures` |
| `analysis/Rad_Dam_Check/Flux_Cal.py` | the full-beam flux for the recorded PIN-diode reading | by hand |

## Code structure

```
.
├── Makefile              the whole build: figures, the three PDFs, make check, make robustness
├── README.md             this file
├── analysis/             data reduction and figure generation; see analysis/README.md
│   ├── environment.yml   the conda environment env-2021-LLDP, matplotlib pinned
│   ├── common/
│   │   ├── acs_style.py  house figure style: column widths, 8 pt Arial, line widths, save_fig()
│   │   └── utils.py      outlier cuts and group averaging, plus beamline-side readers
│   ├── SAXPCS_8id/       the 2021-3 temperature series
│   │   ├── data_2021-3/                  committed reduced data (127 MB): 85 averaged NeXus files
│   │   │                                 and trace_2021-3.csv, one row per acquisition
│   │   ├── temperature_series_2021-3.py  figure script: Fig. 2, Figs. S4-S8
│   │   ├── arrest_vs_temperature.py      figure script: Fig. 3, Fig. S9
│   │   ├── calibration_figure.py         figure script: Fig. S2
│   │   ├── general_model_check.py        robustness numbers, no figure
│   │   ├── contrast_2021-3.py            speckle contrast, no figure
│   │   ├── xpcs_fit.py       shared: the Equation 2 model, its fits, beta from a static standard
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
    --> reduce_babnigg202110.py      groups of about 100 acquisitions, two outlier cuts
analysis/SAXPCS_8id/data_2021-3/               committed; everything below runs from here
    --> temperature_series_2021-3.py  --> Fig. 2, Figs. S4-S8
    --> arrest_vs_temperature.py      --> Fig. 3, Fig. S9
    --> general_model_check.py, contrast_2021-3.py  --> numbers on stdout
abs_xsec.py constants                 --> calibration_figure.py   --> Fig. S2
analysis/Rad_Dam_Check/cluster_results/ --> g2_SAXPCS_Rad_Cali.py --> Fig. S3
manuscript/*.tex and the figures above --> make papers --> manuscript/build/*.pdf
```

The two temperature-series scripts share their building blocks instead of
copying them. `xpcs_fit.py` holds the model and the fit, `series_2021_3.py` the
clock and the SAXS, `nexus_read.py` the file readers, and `abs_xsec.py` the
absolute scale. Every figure script imports `common/acs_style.py`, so all the
figures share one style.

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
mounts only on the beamline GPU hosts (amber, amazonite, adamite).

| Script | Out-of-repo location |
|---|---|
| `analysis/SAXPCS_8id/nexus_convert_babnigg202110.py` | `/gdata/s8id-dmdtn/2021-3/babnigg202110` (the Data Exchange archive) and `..._nexus` (the NeXus tree it writes) |
| `analysis/SAXPCS_8id/restore_dx_metadata_2021-3.py` | the same pair; checks that every file carries its monitors |
| `analysis/SAXPCS_8id/make_qmap_babnigg202110.py` | the same pair, the detector blemish map and the pySimpleMask source under `/home/beams*/8IDIUSER` |
| `analysis/SAXPCS_8id/sweep_stuck_pixels_2021-3.py` | the NeXus tree |
| `analysis/SAXPCS_8id/run_babnigg202110_correlations.sh`, `stop_babnigg202110_correlations.sh` | the NeXus tree; drive `boost_corr_bin` on a GPU host |
| `analysis/SAXPCS_8id/reduce_babnigg202110.py` | `/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results`, the per-acquisition results `data_2021-3/` is averaged from. It also imports `average_ranges.py`, which is not part of this repository. |
| `analysis/common/utils.py`, beamline-side helpers | wherever the caller points them; no default |

All paths are hard-coded at the top of each script. Every other script anchors
its input and output on its own directory and runs from any working directory.

### Where the 8-ID-I data live

| Path | What it holds |
|---|---|
| `/gdata/s8id-dmdtn/2021-3/babnigg202110` | Raw detector frames and the 2021 Data Exchange metadata. The archival copy, never written to. |
| `/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus` | NeXus metadata rewritten for the current `boost_corr_bin`, the q-map, and `reprocess_results/`, the per-acquisition correlation results. The `.bin` files are symlinks into the archive. |
| `/gdata/s8id-dmdtn/2025-1/qzhang202503/flux_control_results` | The four flux-control result files, E0171–E0174 (2025). Their raw frames are in the DM experiment qzhang202503, archived on 23 Aug 2025; see `analysis/Rad_Dam_Check/README.md`. |

8-ID-I moved to NeXus around 2024 and `boost_corr_bin` reads only NeXus
metadata, which is why the 2021 beamtime was converted
(`nexus_convert_babnigg202110.py`, manifests alongside it). The series is
correlated against `babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm2.hdf`,
which masks the two stuck pixels (30, 340) and (63, 166).

## Source data integrity

No data file is modified by producing a figure. Every `h5py.File` call in the
figure scripts and in the modules they import opens its file with mode `'r'`;
the only writing opens in `analysis/` are in `make_qmap_babnigg202110.py`,
which writes a new q-map on beamline storage and is never run by `make`.

Nothing is rescaled between the calibration and the plot. `abs_xsec_coef()` is
built from millimetres, so every absolute intensity axis in this work is
**mm⁻¹**, plotted as the coefficient comes out. (The literature more often
quotes cm⁻¹, which is ten times larger.)

## Submission target

PNAS Research Report, Direct Submission (`PNAS_author_guide.pdf`): preferred
length six pages, about 4,000 words, 50 references and four medium figures;
abstract at most 250 words; significance statement at most 120 words; at least
three keywords; SI Appendix as a single PDF.
