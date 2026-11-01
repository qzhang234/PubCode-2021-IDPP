# Dynamic Arrest Follows Assembly on a Temperature-Scaled Common Clock in an Elastin-like Polypeptide

Everything needed to reproduce the PNAS submission: the reduced experimental
data, the analysis code that turns it into every data-derived figure, and the
LaTeX sources for the manuscript, the SI Appendix and the cover letter.

The measurement is small-angle X-ray photon correlation spectroscopy (SA-XPCS)
at APS beamline **8-ID-I**, October 2021 (beamtime 2021-3): four loadings of
the elastin-like polypeptide (VPAVG)30 held at 28.00, 28.25, 28.50 and
28.75 °C. Correlating the detector frames gives the fast and arrested
fractions at five length scales; averaging the same frames gives I(Q).

The structural half of the original Letter -- the 12-ID-B SAXS/WAXS, the
beta-sheet-like peaks, the thermal-cycling reversibility and the 2022-1 data --
was split out on 28 Oct 2026 into
[PubCode-2021-IDPP-Reversibility](https://github.com/qzhang234/PubCode-2021-IDPP-Reversibility).
The combined version is recoverable from the tag `pre-split`.

```
make env        # create the conda environment, once
make figures    # regenerate every data-derived figure from the reduced data
make papers     # compile main.pdf, si.pdf, cover_letter.pdf into manuscript/build/
make all        # figures + papers
make check      # report undefined citations / references / missing figures
make robustness # recompute the SI 5.3 and 8.2 robustness numbers (no figure)
make clean      # remove LaTeX aux files, keep the PDFs
```

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
| `analysis/SAXPCS_8id/reduce_babnigg202110.py` | `/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results`, the per-acquisition results `data_2021-3/` is averaged from |
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

## Layout

| Path | What it is |
|---|---|
| `Makefile` | The whole build: figures, the three PDFs, `make check`, `make robustness`. |
| `README.md` | This file. |
| `PNAS_author_guide.pdf` | The PNAS submission guidelines, the target for this manuscript. |
| `ACS_author_guide.pdf` | ACS figure guidelines; Appendix 2 is the source of the numbers in `analysis/common/acs_style.py`. |
| `analysis/` | Data reduction and figure generation. See `analysis/README.md`. |
| `manuscript/` | LaTeX sources and compiled PDFs. See `manuscript/README.md`. |
| `tools/` | `check_refs.py` and `pdf_comments.py`, manuscript utilities. |
| `.loglogin` | Beamline login log, carried over from the acquisition account. Not used by anything here. |

## Submission target

PNAS Research Report, Direct Submission (`PNAS_author_guide.pdf`): preferred
length six pages, about 4,000 words, 50 references and four medium figures;
abstract at most 250 words; significance statement at most 120 words; at least
three keywords; SI Appendix as a single PDF.
