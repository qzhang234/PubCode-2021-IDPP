# Figure, data and argument review

**Manuscript:** ELP gelation paper (Nano Letters submission)  
**Date:** 2026-08-23


Companion to `introduction_citation_audit.md`. Covers the user's questions 2 and 3: whether every figure and every datum in it is explained, whether the results support the argument, and whether figure cross-references are correct.


---

## A. Independently re-verified findings

Everything in this section I checked by hand against the deposited data and code, not via an agent. These are the items that would change the paper's claims.


### A1. BLOCKER - the 10 C reversibility comparison uses the same capillary twice

`analysis/SAXS_12id/Read_12ID_SAWAXS.py` lines 39-50 map raw files to the four Figure 2 curves:

| Figure 2 curve | raw file | aliquot |
|---|---|---|
| Reference, 10 C   | `SPA1_10C_00070.avg`  | PA1 |
| Measurement, 10 C | `SPA1_10C_00071.avg`  | **PA1 - the same one** |
| Reference, 30 C   | `SPA1B_30C_00078.avg` | PA1 |
| Measurement, 30 C | `SPA3B_30C_00081.avg` | PA3 |

At 30 C the pair is genuinely two aliquots (PA1 vs PA3). At 10 C it is one aliquot measured twice
(scans 70 and 71). The raw directory shows the naming convention unambiguously - PA1 and PA3 are the two
aliquots, and a trailing B/C denotes a repeat scan of the same one:

```
SPA1_10C_00070   SPA1_10C_00071   SPA1B_30C_00078  SPA1C_30C_00079
SPA3_10C_00073   SPA3B_10C_00074                   SPA3B_30C_00081
```

**The genuine second aliquot at 10 C exists and is not used.** Buffer-subtracting and comparing:

| 10 C pair | q 0.008-0.02 | q 0.02-0.05 | q 0.05-0.2 |
|---|---|---|---|
| PA1/70 vs PA1/71 (**as plotted now**) | 1.6 % | 0.5 % | 0.5 % |
| PA1/70 vs PA3/73 (**the real comparison**) | **44 %** | **10 %** | **10 %** |
| PA3/73 vs PA3/74 (PA3 repeatability) | 18.6 % | 1.1 % | 0.6 % |

So the sentence at `main.tex:93` - "At 10 C, the two profiles overlap ... this profile is recovered after
repeated cycling" - currently rests on measurement repeatability, not on reversibility. With the correct
file the two aliquots differ by ~44 % at low q. This is the evidence for the paper's title claim.

**Required action (authors' call - I have not changed the science):** confirm the aliquot-to-label mapping,
regenerate `Merged_Measurement_10C.csv` from `SPA3_10C_00073` (or `_00074`), redraw Figure 2, and soften the
reversibility wording to match what the corrected comparison shows. Note PA3's own low-q scan-to-scan spread
is 18.6 %, so some of the 44 % is sample//beam variability rather than irreversibility - which is itself
worth reporting.

### A2. MAJOR - each plotted "elapsed time" is a group average, never disclosed

The ten `Average_B0147_*` files each average a contiguous run range:

```
runs     1- 200 = 200 acquisitions      runs  801- 950 = 150 acquisitions
runs   200- 350 = 151                   runs  951-1050 = 100
runs   351- 500 = 150                   runs 1051-1150 = 100
runs   501- 650 = 150                   runs 1151-1250 = 100
runs   651- 800 = 150                   runs 1251-1313 =  63
```

Every "t_w = 5039 s" in the text and figures is therefore an average over 63-200 separate 2 s acquisitions
spanning roughly 6-15 min, labelled with the **start** of that interval. Neither document says so. It
matters physically: averaging g2 over a window in which f is falling can itself broaden the decay and
contribute to the stretched-exponential shape.

**Fix:** state it once in SI Section 2.2 and reference it from the Figure 3, S6, S7, S8 captions; quote
group centres alongside the start times.

### A3. MAJOR - detector clock is internally inconsistent with the data

Text says 52 kHz, 19.2 us per frame, 2 s per acquisition. But 100000 x 19.2 us = 1.92 s, not 2 s. The HDF
metadata records `frame_time = 20.00 us` (50 kHz, exactly 2.000 s), and the multi-tau delay list is built
from that value - so **every tau in the paper sits on a 20 us clock while the text claims 19.2 us**, a 4 %
scale error on all reported relaxation times. 52 kHz is the detector's *maximum* rate (si.tex:87 says so);
the data appear to have been taken at 50 kHz.

**Fix:** reconcile to one clock in main.tex:84, main.tex:141 and si.tex:87.

### A4. MINOR - "the 2 s observation window" overstates the accessible range

The longest multi-tau delay actually plotted is **1.638 s**. "Arrested within the 2 s window" appears in the
abstract and five further places. Replace with "the ~1.6 s maximum delay of the 2 s acquisition".

### A5. MINOR - "five contiguous bins" misdescribes the reduction

The deposited files contain **27** dynamic q bins spanning 0.00376-0.03307 A^-1; the analysis uses the five
lowest (0.00376-0.00827). si.tex:119 reads as though five is all there is. State the selection and its
reason.

### A6. Confirmed correct
- Every number quoted in the text reproduces the deposited fits exactly (f = 0.426-0.716 -> 0.015-0.137;
  p_fast 0.505-0.706; p_slow 0.256-0.490; gamma_fast -2.25/-2.36/-2.11/-1.83/-4.09; 47 min).
- All five plotted fits have reduced chi2 = 0.77-1.89, so the chi2 >= 10 exclusion never fires on shown data.
- Typeset numbering is correct: main Figures 1-3, SI S1-S8, every mention resolves; the Video S1 still is in
  a `center` environment and correctly does not consume an S-number.


---

## B. Are all figures and data explained?

Agent audit, cross-checked against each plotting script.

**Summary:** I read main.tex, si.tex and all six plotting scripts, then rendered every data-derived PDF (plus Setup.pdf, FigureS1 and the Video S1 still) to check what is actually drawn against the captions and text. 26 findings; 4 are blockers.

The most serious: (1) the Figure S8 caption says its curves are the global fits with shared stretching exponents, but g2_grid_SI.py fits each q bin separately with p_fast = p_slow hard-coded to 0.5 — the caption describes a model that was not used; (2) the entire WAXS inset of Figure 2 — the dataset shown (Reference 30 °C only), the red Gaussian+background fits, the red dashed peak markers and their labels — is not mentioned anywhere, even though the beta-sheet spacings come from it; (3) Figure S7 encodes elapsed time by colour in all four panels with no colour key in the figure and no statement in the caption, leaving panels (b) and (c) fully unlabelled; (4) every plotted "elapsed time" is really the start of an average over 54-173 separate 2 s acquisitions spanning the following ~10-15 min (average_ranges.py + /xpcs/average/file_list), which neither the manuscript nor the SI ever says — the text presents them as single times.

Other substantive gaps: Figure 3's caption calls all error bars "fit uncertainties" when panel (b)'s are the measured g_2 uncertainties; Figure S4b never states its q value (it is 0.00304 A^-1, not one of the five main-experiment bins); unexplained grey guide lines (Fig 3c, S7a,d) and the dotted p = 1 line (S7a); the absolute cm^-1 axis of Figure 3a is never tied to the SI calibration (which names only S6); Figure S6's white colour-bar ticks are unexplained and its plasma normalisation differs from Figure 3's, so the same time is a different colour in the two figures; Figure 3's shared key lists a "0 s" dataset that appears in only one of the three panels; and Figure 1's detection schematic (2D pattern, colour bar, q-ring, beamstop mask), its red/blue schematic curves, its cell photographs and its "~10 um" annotation are all undescribed — the last risks being read as a domain size when the XPCS probes 76-167 nm. Figure S1's six numbered callouts have no parts key. Nothing in this list is cosmetic; each is something a reader cannot resolve from the caption, main text or SI.


### BLOCKER (4)


**manuscript/si.tex:251 (Figure S8 caption) vs analysis/SAXPCS_8id/g2_grid_SI.py:42,46-49,89-98,144**  
The Figure S8 caption states "solid curves are global fits with shared stretching exponents at each elapsed time." The curves drawn are not those fits: g2_grid_SI.py fits each q bin independently with scipy curve_fit and with p_fast = p_slow = 0.5 HARD-CODED (module-level `p1, p2, contrast = 0.5, 0.5, 0.135`; only tau_fast, f, tau_slow are free). Nothing is shared and no exponent is fitted, so a reader comparing S8 to Figure 3 is looking at a different model than the one described in Eq. 2 and SI Section 4.

*Fix:* Either regenerate S8 from the global-fit parameters already computed in saxpcs.py (preferred, so main text and SI show the same model), or change the caption to: "Solid curves are two-mode fits to each q bin individually with the stretching exponents fixed at p_fast = p_slow = 0.5 and contrast fixed at 0.135; they are not the global fits of Figure 3."


**manuscript/main.tex:98 (Figure 2 caption) and main.tex:93 (main text) vs analysis/SAXS_12id/Plot_12ID.py:104-161**  
Figure 2 contains a large inset (top right, ~half the panel) that is never mentioned in the caption or anywhere in the main text. Four distinct visual elements in it are unidentified: (i) the black open circles are Reference 30 °C data only, for q > 0.22 A^-1 — the reader is told in the text that BOTH aliquots develop the peaks and cannot tell which one is plotted; (ii) the two red solid curves are Gaussian + linear-background fits over 0.48-1.00 and 1.18-1.70 A^-1; (iii) the two red dashed vertical lines mark the fitted peak centres; (iv) the red "0.72 A^-1" / "1.38 A^-1" labels are the fitted centres, not read-off maxima. The 4.5 and 8.8 A spacings that carry the beta-sheet claim come from these fits.

*Fix:* Add to the Figure 2 caption: "Inset: high-q region of the Reference profile at 30 °C (open circles). Red curves are independent Gaussian-plus-linear-background fits over 0.48-1.00 and 1.18-1.70 A^-1; red dashed lines and labels mark the fitted peak centres, q = 0.72 and 1.38 A^-1 (d = 8.8 and 4.5 A)."


**manuscript/si.tex:243 (Figure S7 caption) vs analysis/SAXPCS_8id/saxpcs.py:572,575,608,654,657**  
In all four panels of Figure S7 the marker/line colour encodes elapsed time, but the figure contains no elapsed-time key of any kind and the caption never says colour means time. The only two legends in the figure key marker SHAPE (p_fast/p_slow in (a), gamma_fast/gamma_slow in (d)); panels (b) and (c) have no legend at all, so five differently coloured tau(q) data sets and five coloured dashed fit lines are completely unlabelled. A reader cannot assign any point or line in S7 to an elapsed time.

*Fix:* Add a shared elapsed-time key to Figure S7 (the same one-row figure legend used in Figure 3 and S8), or at minimum add to the caption: "Colours denote elapsed time as in Figure 3 (5039, 5969, 6595, 7221 and 7863 s from dark to light)."


**analysis/SAXPCS_8id/average_ranges.py:20-31 and data/Average_B0147_*_00801_00950_results.hdf vs manuscript/main.tex:104,111,123 and si.tex:95,164**  
Every plotted "elapsed time" (Figures 3a-c, S6, S7, S8) labels an AVERAGE over 54-173 separate 2 s acquisitions spanning the following ~10-15 minutes, time-stamped with the acquisition time of the FIRST run in the range. The manuscript and SI describe these as single times ("at t_w = 5039 s", "the correlation functions at each elapsed time"), and neither document states anywhere that the plotted I(Q) and g2 are range averages, how many acquisitions each contains, what time window each spans, or that outlier rejection was applied. A reader interprets each g2 as a 2 s snapshot when it is a ~15 min average.

*Fix:* Add to SI Section 4 (and reference it from the Figure 3, S6, S7 and S8 captions): "Each plotted profile and correlation function is the average of 54-173 individual 2 s acquisitions taken at previously unexposed positions over a contiguous run range (after cross-correlation outlier rejection); the quoted elapsed time is the start of that range, which spans roughly 10-15 min." Consider labelling the curves with a range (e.g. "5039-5969 s") instead of a single value.


### MAJOR (11)


**manuscript/main.tex:128 (Figure 3 caption, last sentence) vs analysis/SAXPCS_8id/saxpcs.py:443**  
The caption's blanket statement "Error bars represent fit uncertainties" is wrong for panel (b). The error bars on the g2 points are the measured pointwise g2 uncertainties read from the HDF files, not fit uncertainties; only the panel (c) error bars are fit uncertainties. The same measured error bars appear in Figure S8 with no explanation at all (si.tex:251 mentions no error bars).

*Fix:* Change to: "(b) Error bars are the measured uncertainties of g_2 stored with the multi-tau correlation; (c) error bars are one-standard-deviation fit uncertainties." Add the same sentence about measured g_2 uncertainties to the Figure S8 caption.


**manuscript/si.tex:219 (Figure S4 caption) vs analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py:133**  
Figure S4b plots g2 at a single q bin, but the q value is given nowhere — not on the axes, not in a legend, not in the caption, not in SI Section 5. The script silently takes the first dynamic bin, q = 0.00304 A^-1, which is not one of the five bins of the main experiment (lowest 0.00376 A^-1). Since the whole purpose of the figure is to show that the main-experiment decays are not flux artifacts, the reader cannot check that the control probes a comparable length scale.

*Fix:* State the q in the caption and on the panel: "(b) g_2 at q = 0.00304 A^-1 (2*pi/q = 207 nm), the lowest q bin of the control configuration."


**analysis/SAXPCS_8id/saxpcs.py:500 (Fig 3c), :569-570 and :650-651 (Fig S7a,d) vs manuscript/main.tex:128 and si.tex:243**  
Grey lines connect the points in Figure 3c and in Figure S7a and S7d. They are pure guides to the eye (per-q in 3c, per-parameter in S7), not fits or model curves, and neither caption mentions them. In a figure whose other lines ARE fits (3b solid = global fits; S7b,c dashed = power-law fits), an unlabelled connecting line invites the reader to take it as a fitted trend.

*Fix:* Add "Grey lines connect points at the same q (Fig. 3c) / the same parameter (Fig. S7a,d) and are guides to the eye only." to both captions.


**analysis/SAXPCS_8id/saxpcs.py:578 vs manuscript/si.tex:243 (Figure S7 caption)**  
Figure S7a contains a grey dotted horizontal line at y = 1.0 that is never identified. It marks the simple-exponential limit p = 1, which is exactly the reference the SI text argues about ("Both shared stretching exponents remain below unity"), so it is load-bearing rather than decorative.

*Fix:* Add to the S7 caption: "The dotted horizontal line in (a) marks p = 1, the simple-exponential limit."


**manuscript/figures/Setup.pdf (Figure 1 artwork) vs manuscript/main.tex:89 (caption) and :84 (text)**  
The right half of Figure 1 — the whole detection schematic — is unexplained. The caption describes only the four-stage sequence and says averaging gives I(q) and correlation gives g2. It never mentions: the 2D speckle pattern with its "Photon/Pixel/Frame" colour bar spanning 10^-5 to 10^-2 (is this measured data? from which sample and temperature? per frame?); the white blob at its centre (beamstop/mask); the dotted elliptical ring and the Q vector arrow drawn on it (the azimuthal q annulus that defines a q bin); or the frame stack labelled t, t+dt0, t+2dt0 with I(Q,t).

*Fix:* Extend the caption, e.g.: "Right: a representative 2D speckle pattern (colour bar: photons per pixel per frame; the white central region is the beamstop mask). The dotted ring and Q arrow indicate one azimuthal q bin. Successive frames t, t+dt0, t+2dt0 (dt0 = 19.2 us) are averaged to give I(Q) and correlated to give g2(Q,dt)."


**manuscript/figures/Setup.pdf (Figure 1 artwork, g2 and I(Q) sketch panels) vs manuscript/main.tex:89**  
Both schematic sub-plots in Figure 1 contain a red curve and a blue curve whose meaning is never stated. In the g2(Q,dt) sketch the blue curve decays and the red curve stays high (and is drawn slightly RISING, which reads as increasing correlation with delay); in the I(Q) sketch the red curve lies above the blue. Red and blue are used elsewhere in the same figure for the Heat and Cool arrows and for the frames around the two cell photographs, so the reader has to guess that red = 30 °C (assembled/arrested) and blue = 6-10 °C (dispersed/mobile).

*Fix:* Add: "In the schematic I(Q) and g2 panels, red denotes the high-temperature (30 °C) assembled state and blue the low-temperature dispersed state; the curves are illustrative, not data." Also redraw the red g2 curve as flat rather than rising.


**manuscript/figures/Setup.pdf (two cell photographs, Heat/Cool arrows, "3 mm" and "~10 um" annotations) vs manuscript/main.tex:89**  
The centre-left of Figure 1 shows two photographs of the aluminium cell, one framed in blue and one in red, joined by red "Heat" and blue "Cool" arrows, with a yellow "3 mm" dimension on the red-framed one and a "~10 um" dimension over the magnified stage-iv panel. None is mentioned in the caption. "~10 um" is particularly hazardous: it is presumably the beam footprint / the size of the depicted field of view, but placed directly under the coarsened-morphology cartoon it reads as the domain size, contradicting the 76-167 nm length scales the XPCS actually probes.

*Fix:* Add to the caption: "Photographs show the aluminium cell below (blue frame) and above (red frame) the transition; the 3 mm dimension is the liquid-chamber diameter. The ~10 um bar on stage iv indicates the beam footprint over which the cartoon morphology is sampled, not a domain size."


**analysis/SAXPCS_8id/saxpcs.py:409 and saxs_evolution.py:147 vs manuscript/main.tex:128 (Fig 3a caption) and si.tex:118**  
Figure 3a and Figure S6 plot I(Q) with an ABSOLUTE axis, cm^-1, i.e. the differential scattering cross section per unit sample volume obtained from the flux/transmission calibration of Eqs. S2-S3. The Figure 3 caption says only "Buffer-subtracted SAXS profiles" and never states that the scale is absolute or how it was obtained; the SI sentence that ties the calibration to a figure names only Figure S6 ("The absolute-scale SAXS evolution in Figure S6 used this transmission and flux calibration"), leaving Figure 3a's cm^-1 axis unsourced.

*Fix:* In the Figure 3 caption write "(a) Buffer-subtracted SAXS on an absolute scale, dSigma/dOmega in cm^-1, using the flux and transmission calibration of Figure S3 (SI Section 3.2)", and amend si.tex:118 to read "the absolute-scale SAXS in Figure 3a and Figure S6".


**analysis/SAXPCS_8id/saxs_evolution.py:158 vs manuscript/si.tex:235 (Figure S6 caption)**  
The Figure S6 colour bar carries a set of white horizontal tick marks drawn across it that are not explained anywhere. They mark the elapsed time of each individual averaged acquisition (so the reader can see the times are unevenly spaced and where the Figure 3 subset sits) — information the caption would otherwise not convey at all.

*Fix:* Add to the S6 caption: "White marks on the colour bar give the elapsed time of each measured profile."


**analysis/SAXPCS_8id/saxpcs.py:118-123,320 vs analysis/SAXPCS_8id/saxs_evolution.py:105 and manuscript/si.tex:235**  
Figure 3 and Figure S6 both encode elapsed time with the plasma colormap but use different normalisations, so the same time is a different colour in the two figures: in Figure 3 the five times are mapped to evenly spaced colormap positions (5039 s = dark violet), while in S6 time is mapped linearly over 0-7863 s (5039 s = orange). The S6 caption explicitly invites the cross-figure comparison ("The subset from 5039 to 7863 s corresponds to the SAXS and SA-XPCS data analyzed in Figure 3") without warning that hues do not correspond. The 0 s dataset compounds this: it is black in Figure 3 but dark purple in S6.

*Fix:* Add to the S6 caption: "Colours here are scaled linearly over the full 0-7863 s range and therefore do not match the discrete colours used for the same times in Figure 3; read times from the colour bar."


**analysis/SAXPCS_8id/saxpcs.py:404,539 vs manuscript/main.tex:128 (Figure 3 caption)**  
The shared key along the bottom of Figure 3 lists six elapsed times (0, 5039, 5969, 6595, 7221, 7863 s) and is placed as a figure-level legend, implying it applies to all three panels. The 0 s entry appears only in panel (a); panels (b) and (c) contain only the last five times. A reader will hunt for the black 0 s correlation function in (b) and the black 0 s point in (c) and not find them. The caption also never identifies the black 0 s curve or the blue-square "6 °C ref" curve in panel (a), nor says that (a) shows only the first and last five profiles of the series shown in full in Figure S6.

*Fix:* State in the caption: "(a) shows the 6 °C reference (blue squares), the profile at the start of the hold (0 s, black) and the five later profiles; the correlation data in (b,c) are the five later times only. The full profile series is in Figure S6."


### MINOR (7)


**analysis/SAXS_12id/Plot_12ID.py:95,100 vs manuscript/main.tex:98 (Figure 2 caption) and si.tex:105**  
The black dashed line in the main panel of Figure 2 and its "~Q^-3.2" annotation are never identified in the caption: the reader is not told the line is a fit (rather than a drawn guide slope), that it was fitted to the Reference 30 °C curve only, or over what q range. The range 0.012-0.040 A^-1 exists in the SI but the caption does not point there, and the exponent is the number quoted in the abstract.

*Fix:* Add to the caption: "The dashed line is an unweighted power-law fit to the Reference 30 °C profile over 0.012 < q < 0.040 A^-1, giving I ~ q^-3.2."


**analysis/SAXPCS_8id/saxpcs.py:706-711 vs manuscript/si.tex:211 (Figure S3 caption)**  
Figure S3b cannot be read quantitatively. The y axis is labelled "Number of Photons" (a count, plotted with a 1e9 multiplier) while the caption calls the quantity "Photon flux" (a rate); the x axis is labelled "Upstream Ion Chamber" with no units at all (it is the 1e9-scaled ion-chamber readout of the calibration CSV). The legend's fit equation "6.25e+04 x Up_IC -1.28e+07" therefore has uninterpretable units.

*Fix:* Either relabel the caption to match the axes ("photon counts recorded by the calibrated PIN diode versus the upstream ion-chamber readout; the incident flux follows as F = (CAL_A x Up_IC)/t_exp + CAL_B") or state the units of the ion-chamber readout in the caption.


**analysis/SAXS_12id/Guinier_Plot.py:23-24,72-79 vs manuscript/si.tex:227 (Figure S5 caption) and si.tex:157-161**  
Three things a reader needs for Figure S5 are missing: (i) the error bars on the red points are never explained (they are sigma_I/I, the propagated .avg intensity uncertainty converted to the log scale); (ii) the Guinier fit range is never given numerically — only q_max*R_g = 1.246 is quoted, from which q_max can be back-computed but q_min (0.020 A^-1) cannot; (iii) I_0 = 1.725 +/- 0.023 is quoted with no units while the axis reads ln[I(Q)] (arbs.), so its scale is meaningless to the reader (the script's console output labels it cm^-1, which conflicts with the axis).

*Fix:* Caption: "Weighted linear fit of ln I(q) against q^2 over 0.020 <= q <= 0.055 A^-1 (q_max*R_g = 1.246), giving R_g = 22.65 +/- 0.50 A and I_0 = 1.725 +/- 0.023 in the arbitrary units of the reduced 12-ID-B profiles. Error bars are the propagated intensity uncertainties, sigma_I/I."


**manuscript/figures/FigureS1_Sample_Cells.png vs manuscript/si.tex:182 (Figure S1 caption)**  
The exploded view in Figure S1a carries six circled callout numbers (1-6) and Figure S1b carries two dashed-circle magnification callouts, but no parts key is given in the caption or SI text. The caption names the chamber, the polycarbonate windows, the O-rings and the bolted caps, so the reader can guess some of the numbers but cannot map them, and item 6 (a small cylindrical part) is not described at all.

*Fix:* Add a key to the caption, e.g. "(1) cell body with 3 mm liquid chamber, (2) end caps, (3) bolts, (4) O-rings, (5) 0.127 mm polycarbonate windows, (6) fill plug" (author to confirm the parts), and say what the dashed circles in (b) magnify.


**analysis/SAXPCS_8id/saxpcs.py:532 and :664 ("Elapsed Time (s)") vs manuscript/main.tex:104,121,123 and si.tex:95**  
The x axis of Figure 3c and Figures S7a,d and the S6 colour-bar label all read "Elapsed time", while the entire text discusses the same quantity as the waiting time t_w defined from the RTD reaching 30 °C. No caption states that the plotted elapsed time is t_w. The definitions are not automatically identical: in the scripts the origin is the acquisition start time of the first B0147 averaged file (saxpcs.py:337 `t_ref = start_times[first_file]`), not the RTD event.

*Fix:* Either relabel the axes/colour bar as "Waiting time, t_w (s)" or add "Elapsed time is the waiting time t_w defined in the text (zero when the RTD reached 30 °C, coincident with the first acquisition of the isothermal series)" to the Figure 3, S6 and S7 captions.


**analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py:113,138 vs manuscript/si.tex:219 (Figure S4 caption)**  
Figure S4a plots Q in nm^-1, the only figure in the paper that does not use A^-1, and neither the caption nor SI Section 5 notes the change. Since the figure exists to show that the main-experiment behaviour is not flux-induced, the reader has to convert (0.02-0.6 nm^-1 = 0.002-0.06 A^-1) before they can see that the control covers the same q range as Figure 3a. The y axis "Intensity (photon/pixel/frame)" is also on a different (non-absolute) scale from Figure 3a's cm^-1, which the caption does not mention.

*Fix:* Plot in A^-1 for consistency, or add to the caption: "q is shown in nm^-1 (0.02-0.6 nm^-1 = 0.002-0.06 A^-1, i.e. the q range of Figure 3a); intensities are detector counts per pixel per frame, not on the absolute scale of Figure 3a."


**analysis/SAXPCS_8id/saxs_evolution.py:126,137 vs manuscript/si.tex:235 (Figure S6 caption)**  
The Figure S6 caption does not say the profiles are buffer-subtracted or which buffer measurement was used, although the Figure 3a caption does say "buffer-subtracted" for the same processing. The subtracted buffer (D0138) was recorded ~5 h before the isothermal series and its file name records a different temperature (Buffer_034C), which the reader has no way to learn.

*Fix:* Add "Profiles are buffer-subtracted on the absolute scale, using the same buffer measurement as Figure 3a" to the S6 caption, and state the buffer measurement conditions once in SI Section 3.2.


### NIT (4)


**analysis/SAXPCS_8id/saxpcs.py:559 vs manuscript/si.tex:243**  
Figure S7's panel labels run out of reading order: (a) top-left, (b) top-right, (d) bottom-left, (c) bottom-right. The caption lists them a, b, c, d, so a reader following the caption jumps from the bottom-right panel back to the bottom-left one.

*Fix:* Swap the two bottom panels (put tau_slow bottom-left) or relabel so the letters run left-to-right, top-to-bottom.


**analysis/SAXS_12id/Plot_12ID.py:156 vs :63**  
The Figure 2 inset y axis is labelled "I(Q)" with no units while the main panel is "I(Q) (arbs.)", so the inset's 0.02-0.04 tick values look like an absolute intensity when they are the same arbitrary units.

*Fix:* Label the inset axis "I(Q) (arbs.)" as well, or state in the caption that the inset uses the same arbitrary intensity units.


**analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py:41,84-88**  
The Figure S4 legends list the four fluxes in the non-monotonic order 3.3, 2.7, 4.9, 3.9 x 10^11 ph/s because the entries follow the nominal attenuator order a7, a9, a11, a14 rather than the recovered flux. A reader scanning for a dose trend has to re-sort the legend mentally.

*Fix:* Sort ORDER by the recovered flux so the legend reads 2.7, 3.3, 3.9, 4.9 x 10^11 ph/s.


**manuscript/figures/Setup.pdf (stages i-iv, molecular inset) vs manuscript/main.tex:89**  
Two smaller unexplained encodings in Figure 1: (i) the protein-rich background changes colour across the four stages (pale pink in i-ii, saturated red in iii, dark brick red in iv), which presumably denotes increasing protein concentration in the dense phase but is never stated; (ii) the dashed-box molecular inset (residue colour key plus a "~4 nm" dimension on the VPAVG chain) is never mentioned in the caption, and the ~4 nm chain extent is never related to the measured apparent R_g of 22.65 A.

*Fix:* Add "deepening red indicates increasing protein concentration in the dense phase" and one clause introducing the molecular inset ("lower left: VPAVG repeat with residues coloured as keyed; the ~4 nm bar is the extended length of one repeat unit").


---

## C. Do the results support the argument?

Referee-style assessment. Per your instruction, n=1 and fixed beta are treated as given; the findings below concern only disclosure and framing.

**Summary:** The central claim is largely supported by the data, but the argument-evidence chain has one hard break and several framing gaps that a Nano Letters referee would act on.

What holds up. The dynamic-arrest result is real and robust. Rerunning the deposited global fits reproduces the reported numbers exactly (f = 0.426-0.716 at t_w = 5039 s falling to 0.015-0.137 at 7863 s), all five fits have reduced chi2 of 0.77-1.89, no parameter rails against a bound, and the monotonic decrease in f survives every reasonable perturbation I tried (beta = 0.115/0.155, per-q beta from the measured short-delay contrast, and a 10 q-bin fit). The disappearance of the 8.8 A WAXS feature on cooling is unambiguous (Gaussian amplitude 0.0204 +/- 0.0010 at 30 C versus a negative amplitude at 10 C). The 10 C Guinier state genuinely shows no aggregate upturn. tau_slow is hedged correctly and consistently in four separate places and is not over-interpreted anywhere; the SI's statement that 1-f is neither a volume fraction nor the infinite-time nonergodicity factor is exactly right. SI Section 2.3 on replication is a model of honest disclosure.

The break. The reversibility claim — the paper's title claim — rests on Figure 2's 10 C overlap, and in the deposited code both 10 C traces come from the same capillary (SPA1_10C_00070 and _00071, differing by 0.41%), while the 30 C traces come from two different capillaries. The genuine second aliquot at 10 C exists in the raw directory and differs by 9-10% in the Guinier region with 1.7-6.4x more low-q intensity. Until the file-to-aliquot mapping is resolved and stated, the 10 C overlap is measurement repeatability, not cycled-versus-fresh agreement. Compounding this, all 12-ID 10 C scans precede all 30 C scans, so no post-heating cool-down was measured at all in that beamtime.

The framing gaps. Given that n=1 and the two-experiment design are unavoidable, the question is disclosure, and the disclosure is uneven: both caveats are handled well in the SI and in one main-text limitations paragraph, but the abstract, the framing paragraph and the summary all speak as if one system, and the "Three limitations" paragraph omits replication entirely. The abstract also escalates to "amyloid-like" four times on the strength of two broad rings whose fitted widths correspond to correlation lengths of only ~19 and ~25 A. Fixing beta is physically correct, but the manuscript never tells the reader why, never says how 0.135 was obtained, and never shows the sensitivity — which matters because the specific endpoint values (0.01-0.14) shift to 0.08-0.20 at beta = 0.155.

Also actionable: two undisclosed selections (5 of 10 elapsed times, 5 of 27 deposited q bins — and the discarded higher-q bins fit cleanly and show arrest is markedly less complete at 45 nm, which is a result worth having); the fact that each "elapsed time" is really a 7-16 minute group average labelled by its start; a chi2 < 10 criterion that is disclosed without numbers and is inert on the analyzed data while the real exclusion is justified on different grounds; and a Figure S8 caption that describes global shared-exponent fits when the code plots per-q fits with the exponents nailed to 0.5.

Nearly all of this is fixable with text, numbers already in hand, and one figure regeneration. Only the aliquot-provenance issue requires the authors to go back to their records.


### BLOCKER (3)


**analysis/SAXS_12id/Read_12ID_SAWAXS.py lines 39-42 and 47-50; main.tex line 93; si.tex line 90; Figure 2**  
The reversibility claim rests on overlap between a cycled "Measurement" aliquot and a refrigerated "Reference" aliquot at 10 C, but in the deposited code the two 10 C curves come from the SAME capillary. Reference_10C = SPA1_10C_00070 and Measurement_10C = SPA1_10C_00071 (same sample label PA1, consecutive scan numbers), whereas the 30 C pair is Reference_30C = SPA1B_30C_00078 and Measurement_30C = SPA3B_30C_00081 (different capillaries, PA1 vs PA3). The 10 C "overlap" is therefore scan-to-scan repeatability of one capillary, not cycled-vs-fresh agreement, and the 10 C and 30 C "Measurement" traces are not the same physical sample.

*Fix:* Resolve the file-to-aliquot mapping and state it explicitly. If PA1 and PA3 are the two aliquots, Measurement_10C must be regenerated from SPA3_10C_00073 or SPA3B_10C_00074 and Figure 2 redrawn; the low-q excess of that aliquot at 10 C must then be shown and discussed (it may indicate incomplete mesoscale recovery on cooling, which would require softening "the low-temperature profile is recovered" to "the molecular-scale profile is recovered; the lowest-q intensity of the cycled aliquot remains X% above the reference"). Add to SI Section 2.1 a sentence naming the file/scan for each of the four traces, e.g. "Reference = capillary PA1 (scans 00070, 00078); Measurement = capillary PA3 (scans 00074, 00081)."


**main.tex lines 72-73 (abstract) and 136 (summary)**  
The abstract asserts "coexistence of thermally reversible amyloid-like local order with progressive mesoscale arrest demonstrates..." and the summary says "(VPAVG)30 produces a hierarchical protein-rich assembly with reversible beta-sheet-like local order and progressively suppressed mesoscale fluctuations" — both present one system. In fact the structural and dynamical data come from different samples, different beamlines, and beamtimes ~11 months apart, and the arrested 8-ID sample was never shown to carry the WAXS signature at all (8-ID covers only q = 0.0038-0.033 A^-1). Neither the different-sample caveat nor n=1 appears anywhere in the abstract.

*Fix:* Add one clause to the abstract, e.g. "...in a separate isothermal experiment on an independently prepared sample, SA-XPCS at 30 C shows..." and replace "demonstrates" with "indicates". In the summary (line 136) replace "produces a hierarchical protein-rich assembly with reversible beta-sheet-like local order and progressively suppressed mesoscale fluctuations" with "...produces a hierarchical protein-rich assembly that shows reversible beta-sheet-like local order in temperature-dependent SAXS/WAXS and, in a separate isothermal SA-XPCS experiment, progressively suppressed mesoscale fluctuations." Also expand limitation 2 (line 134) to say explicitly that the arrested sample was not itself measured at WAXS q.


**main.tex line 134 ("Three limitations"); si.tex Section 2.3 (line 98); main.tex line 141**  
n = 1 for the isothermal series is disclosed only in SI Section 2.3 and in one sentence at the end of the Experimental Section. It is absent from the "Three limitations" paragraph, which is where a referee looks for the honest inventory, and absent from the abstract. As written, a reader who skims the limitations paragraph comes away believing the only caveats are spectroscopy, separate experiments, and the finite time window.

*Fix:* Promote it. Make it a fourth limitation in main.tex line 134: "Fourth, the isothermal series was acquired from a single loaded sample; the ~1300 fresh-position acquisitions provide spatial sampling within one preparation, and the quoted uncertainties are fit uncertainties, not between-sample variability. The reported trend in f is therefore a single-preparation observation." This costs two lines and removes the strongest referee objection to the framing rather than the data.


### MAJOR (6)


**si.tex line 134; analysis/SAXPCS_8id/saxpcs.py line 128**  
beta = 0.135 is stated as fixed "as determined from the measured short-delay contrast" with no physical justification, no statement of which dataset/delay range determined it, and no sensitivity analysis. Fixing beta is physically correct (it is an instrument parameter), but the manuscript never tells the reader that, and the numbers most quoted in the abstract and conclusions are the ones most sensitive to the choice.

*Fix:* Add to si.tex Section 4: "beta is the instrumental speckle contrast set by the beam coherence and detector pixel size. It is a property of the beamline configuration, not of the sample state, and was therefore held fixed at the value measured from the shortest-delay g2, beta = 0.135 (mean over the five analyzed acquisitions and five q bins at delays < X s). Allowing beta to float per elapsed time drives it to unphysical values above the measured short-delay contrast while improving reduced chi2 by less than 0.15, confirming that it is absorbing model misfit rather than a real change. Repeating the analysis with beta fixed at 0.115 and 0.155 changes f at t_w = 7863 s to 0.00-0.09 and 0.08-0.20 respectively; the monotonic decrease of f at every q is unchanged, so the trend is robust while the absolute endpoint values carry this systematic." Then soften "near-complete shift" (main.tex line 136), which is the one claim that does not survive the beta systematic.


**analysis/SAXPCS_8id/saxpcs.py lines 110-112 (N_LAST = 5, fit_q_indices = [0,1,2,3,4]); si.tex line 125; main.tex line 123**  
Two data selections are made silently. (i) Only the last 5 of 10 averaged isothermal files are fitted. (ii) Only 5 of the 27 q bins present in the deposited reduction are analyzed. si.tex line 125 describes the reduction as if five bins were all that exist ("Azimuthal averaging was performed in five contiguous bins"), which is not what the deposited files contain. The 76-167 nm range quoted in the abstract is thus a chosen window, not the accessible one.

*Fix:* State both selections and, ideally, use the extra information. Reword si.tex line 125 to "The reduction produced 27 azimuthally averaged q bins spanning 0.0038-0.0331 A^-1. Two-step correlation functions with adequate signal at all analyzed elapsed times were obtained for the five lowest-q bins, which are reported here." And consider adding the higher-q bins at the late times to Figure S7 — the fact that f at 7863 s is 0.015 at 167 nm but 0.24 at 45 nm is a real, physically interesting result (arrest is progressively less complete toward smaller scales) that the current window discards, and reporting it strengthens rather than weakens the paper.


**main.tex lines 104, 111, 123, Figure 3 caption; si.tex lines 95, 164**  
Each plotted "elapsed time" is not a time point: it is an average over 63-200 separate 2 s acquisitions spanning 6.6-15.5 minutes, and the label is the START of that interval, not its center. This is never disclosed. It matters physically: averaging g2 over a window in which f is falling steeply will itself broaden the decay and can contribute to the apparent stretched-exponential shape and to the two-step form that the whole analysis rests on.

*Fix:* Add to si.tex Section 2.2: "Each reported elapsed time is a group average over 63-200 consecutive 2 s acquisitions at fresh positions, spanning 6.6-15.5 min; the quoted t_w is the start of the group (group centers are 5504, 6282, 6908, 7542 and 8061 s). Because f varies by less than X within a group, the averaging does not materially broaden the correlation functions." Either quote group centers or give the interval (e.g. "5039-5969 s") in the Figure 3 key and caption, and correct "for up to 7863 s" to "for approximately 8300 s".


**si.tex line 136 ("Global fits with reduced chi^2 >= 10 were excluded"); analysis/SAXPCS_8id/saxpcs.py line 114; si.tex line 164; main.tex line 104**  
The chi2 < 10 exclusion is disclosed as a bare threshold with no numbers, which reads like undisclosed data rejection, and it is in fact inert on the analyzed data. Separately, the real exclusion — the four earlier elapsed times — is justified in the text on a different ground ("insufficient coherent intensity") than the one implemented in the code.

*Fix:* Replace si.tex line 136's bare threshold with the numbers and the correct reason: "Fits with reduced chi^2 >= 10 were excluded. All five analyzed elapsed times satisfied this criterion (reduced chi^2 = 0.77-1.89). At t_w <= 4088 s the correlation functions decayed fully within the accessible delay range with no resolvable second step (g2(tau = 1 s) - 1 < 0.01) and the two-mode model failed (reduced chi^2 = 48-1859); those times are therefore excluded from the fitted trends and appear only in the SAXS evolution of Figure S6." This is defensible and strengthens the arrest narrative — the excluded times are ergodic, not merely noisy.


**si.tex line 251 (Figure S8 caption); analysis/SAXPCS_8id/g2_grid_SI.py lines 42-49 and 94**  
The Figure S8 caption states "solid curves are global fits with shared stretching exponents at each elapsed time." They are not. g2_grid_SI.py performs independent per-q scipy curve_fit calls with p_fast and p_slow HARD-CODED to 0.5, and no sharing across q. So the fits displayed at four of the five q values are not the fits from which the reported f values come.

*Fix:* Regenerate Figure S8 from the saved global-fit results in saxpcs.py (preferred, so the SI shows the fits that produced the reported f), or else correct the caption to "solid curves are independent per-q fits of Eq. S4 with p_fast = p_slow = 0.5 held fixed, shown to illustrate the data quality; the fitted parameters reported in Figure 3c and Figure S7 come from the global shared-exponent fits described in Section 4."


**main.tex lines 73 (abstract, x3), 76, 136; caption of Figure 2 (line 98)**  
The abstract and conclusion escalate from "beta-sheet-like" (which the body carefully justifies) to "amyloid-like", four times, on the basis of two broad WAXS rings and no spectroscopy, no ThT, no CD/FTIR, and no fibril morphology. The peak widths — which are the strongest argument that this is short-range local order rather than amyloid — are never reported.

*Fix:* Report the widths — they help the argument. Add to Figure 2's caption or si.tex Section 3.1: "The two features have FWHM of 0.34 and 0.26 A^-1, corresponding to correlation lengths of only ~19 and ~25 A, i.e. local order extending over a few repeat spacings rather than extended cross-beta crystallites." In main.tex line 102 add one sentence naming the alternative: "A ~4.5 A ring is also produced by generic dense-polypeptide packing, so the assignment rests on the concurrent appearance of both spacings and on their thermal reversibility rather than on either peak alone." Then restrict "amyloid-like" to the introduction's framing and use "beta-sheet-like" in the abstract, the Figure 2 caption and the summary.


### MINOR (7)


**main.tex line 93; analysis/SAXS_12id/Plot_12ID.py lines 106-138**  
"At 30 C, both aliquots reproducibly develop an approximate I(q) ~ q^-3.2 response at low q and peaks at q = 1.38 and 0.72 A^-1" — but only the Reference aliquot is fitted for either quantity, and only the Reference appears in the WAXS inset. Fitting the Measurement aliquot gives noticeably different peak parameters.

*Fix:* Either report both, e.g. "peaks at q = 1.38 and 0.72 A^-1 (Reference) and 1.36 and 0.74 A^-1 (Measurement)", or reword to "both aliquots develop an approximate q^-3.2 low-q response (slopes -3.17 and -3.15) and the same pair of high-q features, at q = 1.38 and 0.72 A^-1 in the Reference aliquot." Add the Measurement 30 C trace to the inset so a reader can see the comparison the text asserts, and note in si.tex Section 3.1 that peak amplitudes are not on a common absolute scale because the WAXS-to-SAXS overlap factor is fitted per profile and its uncertainty is not propagated.


**analysis/SAXS_12id/Read_12ID_SAWAXS.py lines 76-79; main.tex line 93; si.tex line 103**  
Undocumented q-range crops make the 10 C and 30 C profiles cover different q ranges, and this affects two claims. The 10 C data are truncated at q = 1.41 A^-1 (30 C data reach 1.61), so the disappearance of the 4.5 A peak on cooling is documented over only the low-q flank of that peak, and the manuscript's own fit window for that peak (1.18-1.70 A^-1) does not exist in the 10 C data. At low q, the 10 C data start at 0.0081 A^-1 while the 30 C data start at 0.0046 A^-1, so "no detectable mesoscale assembly" at 10 C is established only for length scales below ~78 nm — which barely overlaps the 76-167 nm XPCS window.

*Fix:* Document the crops in si.tex Section 3.1 with the reason for each, and state the accessible ranges: "the 10 C profiles are reported over 0.0081-1.41 A^-1 and the 30 C profiles over 0.0046-1.61 A^-1." Then hedge accordingly in main.tex line 93/102: "the absence of mesoscale assembly at 10 C is established for q >= 0.008 A^-1 (structures below ~78 nm)" and note that the 4.5 A feature's disappearance is documented over the accessible 1.18-1.41 A^-1 interval while the 8.8 A feature is fully resolved at both temperatures.


**main.tex line 84 and line 141; si.tex line 93**  
Internally inconsistent detector timing, contradicted by the deposited data. The text says 52 kHz, 19.2 us per frame, and 2 s per acquisition — but 100,000 x 19.2 us = 1.92 s, not 2 s. The HDF files record frame_time = 20.0 us (i.e. 50 kHz, 2.000 s), and the multi-tau delay list is built from that value, so every tau in the paper is on a 20 us clock while the text says 19.2 us.

*Fix:* Reconcile to one clock. If the acquisition really ran at 50 kHz / 20 us / 2.000 s, correct "52 kHz" and "19.2 us" in main.tex lines 84 and 141 and si.tex line 93. If 52 kHz is correct, state the acquisition as 1.92 s and explain the 20 us value in the metadata. Either way it is a 4% scale on every reported tau.


**main.tex lines 73, 111, 121, 123, 134; si.tex lines 138, 169**  
"The 2 s observation window" is used throughout as the arrest criterion, but the longest measured delay is 1.64 s. The claim "arrested within the 2 s window" overstates the accessible range by ~20%.

*Fix:* Replace "the 2 s observation window" with "the ~1.6 s maximum delay of the 2 s acquisition" (or "the accessible delay range, 2 x 10^-5 to 1.6 s") in the abstract and at each occurrence. This costs nothing and removes a discrepancy any referee can check against Figure 3b.


**main.tex line 104 and line 141; si.tex line 95**  
The stated definition of t_w ("zero when the calibrated RTD readout reached 30 C") does not match the time origin implemented in the analysis, which is the start time of the first isothermal acquisition. On the stated ramp protocol these differ by several minutes, so every quoted t_w may be systematically offset.

*Fix:* Either reconcile the two (if the first acquisition began when the RTD reached 30 C, say so and explain the 21 min gap after the 6 C reference), or state the offset explicitly in si.tex Section 2.2: "t_w = 0 is taken as the start of the first isothermal acquisition, which followed the RTD reaching 30 C by approximately N s; all quoted t_w therefore underestimate the time at 30 C by that constant offset." A constant offset changes no conclusion but is a reproducibility problem as written.


**main.tex line 123; si.tex line 169; Figure S7 panel (d)**  
"the final-time exponent becomes poorly constrained as f approaches zero" mis-describes what Figure S7d actually shows. The final-time gamma_fast is not imprecise — it is precisely a different value, and a referee comparing text to figure will see error bars comparable to the other points.

*Fix:* State the value and the reason plainly: "gamma_fast is close to -2 over most of the analyzed interval (-1.8 to -2.4) but departs sharply at the final time (-4.1 +/- 0.6). At that time f < 0.14 at every q, so the fast mode contributes little amplitude and tau_fast(q) is scattered and non-monotonic; we therefore do not interpret the final-time exponent." Same correction in si.tex line 169.


**si.tex Section 2.1 (line 90) and Section 5 (line 143)**  
Two pieces of provenance a referee needs are missing: (a) the equilibration time at 30 C before the 12-ID exposures, without which the structural state cannot be placed anywhere on the 8-ID isothermal kinetic axis (where arrest took >5000 s to develop); (b) the temperature and physical state of the 20 mg/mL radiation-control sample, without which the flux control cannot be judged as bearing on the arrested state.

*Fix:* Add the 30 C hold time before each 12-ID exposure to si.tex Section 2.1, and add the control temperature/state to Section 5, e.g. "...at 30 C in the assembled state" or "...at 6 C in the dispersed state; the control therefore bounds gross flux artifacts in the [state] but does not test dose-accelerated arrest in the assembled state." Also note the control samples one instant rather than a 47 min evolution — the fresh-position protocol, not this control, is what rules out cumulative dose.


### NIT (2)


**si.tex Section 8.2 (lines 157-161); analysis/SAXS_12id/Guinier_Plot.py lines 23-24**  
The Guinier fit range is not stated in the SI, only q_max*Rg = 1.246. Because usable data exist well below the fit window, a referee will want to know whether the reported Rg depends on where the window starts — this is the standard check for aggregate contamination in a concentrated sample.

*Fix:* State the window and the check: "The Guinier fit used 0.020 <= q <= 0.055 A^-1 (q_max Rg = 1.246). Data at lower q (down to 0.0081 A^-1) lie within 2% of the Guinier extrapolation with no upturn, confirming the absence of detectable aggregates over the accessible range." That converts an unstated choice into positive evidence for the claim.


**main.tex line 82 and line 132**  
Two framing sentences do the heavy lifting for the central claim and both quietly merge the two experiments. Line 82 says "Separately, isothermal SA-XPCS..." — good instinct, but "separately" is ambiguous between "in a separate measurement on the same sample" and "in a separate experiment on a different sample." Line 132 ("the high-temperature state contains thermally reversible beta-sheet-like contacts while the low-q assembly progressively loses mobility") reads as a single state described by two probes.

*Fix:* Line 82: "In a separate experiment on an independently prepared sample, isothermal SA-XPCS at 30 C reveals..." Line 132: "The present measurements add a local structural element from a separate preparation: the high-temperature state of (VPAVG)30 contains thermally reversible beta-sheet-like contacts, and, in the isothermal experiment, the low-q assembly progressively loses mobility." With those two edits plus the abstract clause, the separate-experiment caveat is visible everywhere the coexistence claim is made rather than only in the limitations paragraph.


---

## D. Figure cross-referencing

Every figure is cited, and all hard-coded numbers currently match the float order LaTeX assigns.

**Summary:** Audited all figure/video cross-references in /home/8-id-i/QZ/Papers/PubCode-2021-IDPP/manuscript/main.tex and si.tex against the compiled build/main.aux and build/si.aux, the rendered figure PDFs (panel-letter bounding boxes extracted with pdftotext -bbox), the two photographic PNGs, and the seven generating scripts under /analysis. Twelve findings; three are major.\n\nWhat is clean. (1) Coverage and ordering: every one of the 11 figures plus Video S1 is cited at least once, and every float is cited before it appears — verified in the compiled output, where main.pdf places the Figure 1/2/3 captions at text lines 141/235/428 against first citations at 69/146/258, and every SI float uses [H] after a \\clearpage so source order is pinned. (2) Hard-coded numbering: all twelve hard-coded \"Figure~N\"/\"Figure~S#\" strings match what LaTeX actually assigns — build/main.aux gives fig:setup->1, fig:static->2, fig:xpcs->3, and build/si.aux gives fig:S1->S1 through fig:S8->S8. (3) The Video S1 still is confirmed to be inside a `center` environment, not a `figure` (si.tex:173-175), so it does not step the figure counter; pdftotext of build/si.pdf shows it rendering as \"Video S1. ...\" with no \"Figure S#:\" prefix, and Figures S1-S8 are consequently numbered correctly with no off-by-one. (4) No reference anywhere points to a panel letter that does not exist: Figure 3a/3b/3c and Figure S3a/S3b all resolve correctly, and the (a)-(d) letters on the Video S1 still and the (a)/(b) on Figure S1 match their captions in the actual images.\n\nThe three major findings. First, SI figures are not numbered in order of first mention — Figure S2 is cited at si.tex:79, before Figure S1 at si.tex:93, but the floats are placed S1-then-S2; the two should be swapped. Second, the Figure S8 caption (si.tex:251) states twice that the curves are \"global fits with shared stretching exponents,\" but g2_grid_SI.py fits each q bin independently with p_fast and p_slow hard-coded to 0.5 (line 42) and only three free parameters (line 52) — this also contradicts the procedure the SI describes at si.tex:134, and Figure 3b genuinely does use the global fit, so the two figures are fitted differently. Third, Figure S7's panel letters are confirmed out of reading order: saxpcs.py:557 lays out axp/axf/axg/axs as TL/TR/BL/BR but line 559 calls label_panels((axp, axf, axs, axg)), putting (c) at bottom-right and (d) at bottom-left (verified by PDF bounding boxes). The caption at si.tex:243 names each letter correctly, so there is no factual error — but a reader following (a),(b),(c),(d) zig-zags, and the fix requires reordering the caption as well as the label call, including changing \"panels (b,c)\" to \"panels (b,d)\".\n\nRemaining items are smaller: no \\ref is used anywhere despite 11 defined labels (currently all correct but unprotected); the main-text SI paragraph (main.tex:144) omits SI Sections 6 and 7, which the SI's own Contents list does include; four unqualified \"Figure 2\"/\"Figure 3\" references inside the standalone SI; two Contents entries that bundle Section 2 material under Sections 4 and 5; an undescribed inset in Figure 2; an undefined panel-letter set in the Figure S8 caption; and a nm^-1 abscissa in Figure S4a where everything else uses Å^-1.


### MAJOR (3)


**manuscript/si.tex:79 (first mention of Figure S2) vs si.tex:93 (first mention of Figure S1); floats at si.tex:179 (S1) and si.tex:187 (S2)**  
SI figures are not numbered in order of first mention. Figure S2 (oligonucleotide/OERCA template sequences) is cited first, at si.tex:79 in Section 1.2 "Gene construction"; Figure S1 (sample cells) is not cited until si.tex:93 in Section 2.2. But the floats are placed S1-then-S2, so LaTeX assigns the sample-cell figure S1 and the sequence figure S2 — the reverse of citation order. This is the only ordering violation in either document; every other figure (S3–S8, main 1–3) is cited before it appears and in ascending order.

*Fix:* Swap the two floats: move the oligonucleotide/OERCA block (si.tex:187–205) ahead of the sample-cell block (si.tex:179–184), swap the \label names (fig:S1 <-> fig:S2), and swap the two hard-coded strings — si.tex:79 becomes `(Figure~S1)` and si.tex:93 becomes `(Figure~S2)`. Also update the file name FigureS1_Sample_Cells.png -> FigureS2_Sample_Cells.png if you keep names aligned with numbers. Alternatively, if you prefer to keep the current numbering, move the first mention of the sample cell earlier — but the float swap is the cleaner fix.


**manuscript/si.tex:251 (Figure S8 caption); generating script analysis/SAXPCS_8id/g2_grid_SI.py:42, :46, :52, :89, :141**  
The Figure S8 caption describes fits that the figure does not contain. It says, twice, that the curves are global fits with shared stretching exponents ("Correlation functions and global two-mode fits..." and "solid curves are global fits with shared stretching exponents at each elapsed time"). In fact g2_grid_SI.py fits every q bin independently, with both stretching exponents hard-coded to 0.5 rather than fitted or shared. This also contradicts si.tex:134, which tells the reader that all five q bins are fitted jointly with p_fast and p_slow shared.

*Fix:* Preferred: regenerate FigureS8_g2_Grid.pdf using saxpcs.py's `fit_g2_global` so S8 matches Figure 3b and the procedure described at si.tex:134, and note the change per the CLAUDE.md rule on fitting-model changes. Minimum fix if the figure is frozen: rewrite si.tex:251 to "Correlation functions and two-mode fits at the four additional q values... solid curves are fits performed independently at each q with the stretching exponents fixed at p_fast = p_slow = 0.5," and add a sentence at si.tex:138 stating that Figure S8 uses this fixed-exponent per-q fit rather than the global fit of Section 4.


**analysis/SAXPCS_8id/saxpcs.py:557–559 (figure generation); manuscript/si.tex:243 (Figure S7 caption)**  
Figure S7's panel letters run out of reading order: (a) top-left, (b) top-right, (d) bottom-LEFT, (c) bottom-RIGHT. The caption enumerates (a), (b), (c), (d) in that sequence, so a reader following the caption goes top-left -> top-right -> bottom-right -> bottom-left. Each letter does name the correct panel, so there is no factual mismatch, but the layout is a legibility/production defect that reviewers and ACS production routinely flag.

*Fix:* Change saxpcs.py:559 to `label_panels((axp, axf, axg, axs))` so the letters run in reading order, regenerate FigureS7_Fit_Parameters.pdf, and reorder si.tex:243 to match the new assignment: "(a) Stretching exponents... (b) Effective fast relaxation time versus q. (c) Fitted scaling exponents gamma_fast and gamma_slow. (d) Effective slow relaxation time versus q. Dashed lines in panels (b,d) are power-law fits tau ∝ q^gamma." Note the current "Dashed lines in panels (b,c)" must become "(b,d)" — the two power-law panels are the right-hand column.


### MINOR (5)


**manuscript/main.tex:90, :99, :129 and manuscript/si.tex:183, :204, :212, :220, :228, :236, :244, :252**  
All eleven figure \label{} anchors are defined and never \ref{}'d — every figure number in prose and in captions is a hard-coded string. I verified that all of them are currently correct, but nothing in the build would catch it if a float were inserted, moved, or promoted from a center environment.

*Fix:* Replace the hard-coded strings with \ref: `Figure~\ref{fig:setup}`, `Figure~\ref{fig:S5}`, etc. Panel letters can stay literal (`Figure~\ref{fig:xpcs}a`). If you prefer to keep literal numbers for a journal-submission reason, at minimum add a `make check` grep that compares each hard-coded `Figure~S?N` against the corresponding \newlabel value in build/*.aux.


**manuscript/si.tex:172–176 (Supporting Video block)**  
CONFIRMED as described in the task: the Video S1 still is inside a `center` environment, not a `figure` environment, and its label is hand-typed bold body text rather than a \caption. The implication for numbering is benign — it does not step the figure counter, so the eight real floats correctly number S1–S8 and no SI figure number is off by one. But the block is fragile and stylistically inconsistent, and the risk is amplified by the zero-\ref finding above.

*Fix:* Keep it as a non-float (that is the correct choice — it must not consume an S-number), but make it robust and consistent: wrap the still in `\begin{figure}[H] ... \end{figure}` only if you also add `\addtocounter{figure}{-1}`, or simpler, leave the center environment and wrap the label line in `{\small ...}` to match caption sizing. Add a source comment at si.tex:172 stating explicitly that this block must never become a figure environment because Figures S1–S8 are hard-coded.


**manuscript/main.tex:144 (Supporting Information paragraph) vs manuscript/si.tex:56–67 (Supporting Information Contents) and the SI section structure**  
The main-text Supporting Information paragraph is missing two items that the SI actually contains and that the SI's own Contents list advertises: Section 6 "Code and Reduced-Data Availability" (si.tex:145) and Section 7 "AI-Assisted Analysis-Code Development" (si.tex:148). It also never mentions the sample-cell figure or the gene-construction/oligonucleotide sequences, both of which are SI figures.

*Fix:* Extend main.tex:144 to: "...correlation functions at all measured q values; code and reduced-data availability; AI-assisted analysis-code development; and a linked video of the reversible thermal transition." Optionally add "sample-cell drawings and gene-construction oligonucleotide sequences" after "sample environments" so Figures S1 and S2 are accounted for.


**manuscript/si.tex:105, :164, :235, :251**  
Four places in the SI refer to "Figure~2" and "Figure~3" — main-text figures — without saying so. Inside a standalone SI whose own figures are S1–S8, an unqualified "Figure 3" is ambiguous, and the SI is distributed as a separate PDF.

*Fix:* Write "Figure 2 of the main text" and "Figure 3 of the main text" at all four locations (or "main-text Figure 3" in the captions to keep them short).


**manuscript/si.tex:61 and si.tex:62 (Supporting Information Contents items 4 and 5)**  
Two Contents entries bundle material that does not live in the section they appear to name, so the list reads as a table of contents but does not match the section structure.

*Fix:* Either move the temperature-protocol and replication phrases into the item that names Section 2 ("Sample environments, temperature protocols, X-ray measurement protocols, and experimental replication") and trim items 4 and 5 to "SA-XPCS correlation analysis and global two-mode fitting" and "Radiation-damage mitigation and flux controls"; or add a leading note that the list is thematic rather than section-by-section.


### NIT (4)


**manuscript/si.tex:251 (Figure S8 caption)**  
Figure S8 prints panel letters (a)–(d) in the artwork, but the caption never uses them — it just lists four q values in order and leaves the reader to match them to panels. Most journals require every panel letter shown in the graphic to be defined in the caption.

*Fix:* Rewrite the second sentence as "Panels show (a) q = 0.00489, (b) 0.00602, (c) 0.00714, and (d) 0.00827 Å^-1."


**manuscript/main.tex:98 (Figure 2 caption) and main.tex:93**  
Figure 2 contains an inset that neither the caption nor the text identifies as an inset. The WAXS peak region is drawn in a separate inset axes, but the caption describes the peaks as if they were in the main axes, so a reader is not told where to look.

*Fix:* In main.tex:98, change to "...and peaks at q=1.38 and 0.72 Å^-1 (inset), corresponding to real-space spacings of approximately 4.5 and 8.8 Å", and add a closing sentence "Inset: the WAXS region of the Reference profile at 30 °C with Gaussian-plus-linear-background fits (red)."


**manuscript/si.tex:251 (Figure S8 caption)**  
"Colors denote elapsed times used in Figure~3" is slightly imprecise: Figure 3 keys six elapsed times, Figure S8 shows only five.

*Fix:* Change to "Colors denote the same elapsed times as in Figure 3b,c of the main text." Also update the stale g2_grid_SI.py docstring at lines 3–5 to match line 112.


**analysis/Rad_Dam_Check/g2_SAXPCS_Rad_Cali.py:138; figure referenced at manuscript/si.tex:219 and main.tex:123**  
Figure S4 panel (a) plots Q in nm^-1, while every other figure in the paper and all body text use Å^-1. The caption does not state units, so there is no contradiction, but a reader comparing Figure S4a with Figure 3a or Figure S6 will silently mis-scale by a factor of 10.

*Fix:* Change g2_SAXPCS_Rad_Cali.py:138 to `axs[0].set_xlabel(r'$Q$ ($\AA^{-1}$)')` and convert the plotted abscissa accordingly, then regenerate FigureS4_Flux_Control.pdf. If the figure is frozen, add "Note the nm^-1 abscissa in panel (a)" to si.tex:219.
