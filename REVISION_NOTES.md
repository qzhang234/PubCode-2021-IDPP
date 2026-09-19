# Revision notes

Changes made to the manuscript, SI, cover letter and analysis code during the
pre-submission review pass. Every number quoted below was re-derived from the
committed data by re-running the scripts in `analysis/`; nothing was copied
from the previous draft.

Build state: `make papers` is clean — no undefined citations, references or
figures. Abstract 129 words (limit 150). Main text 2498 words excluding
captions, 3080 including them.

---

## 1. Corrections of fact

### 1.1 Adjacent-$Q$ significance of the fast fraction

This was wrong in the submitted draft.

- **Was:** "adjacent-$q$ differences exceed $2\sigma$ at $t_w=5039$ s and $8\sigma$ or more thereafter"
- **Now:** "adjacent-$Q$ differences of 1.8--2.9$\sigma$ at $t_w=5039$ s and $8\sigma$ or more thereafter; across the full $Q$ range the increase is $9\sigma$ or more at every analyzed time"

Refitting the five groups from the committed HDF files gives, at
$t_w = 5039$ s, adjacent differences of 2.0, 2.2, 2.9 and **1.8**$\sigma$. The
last pair fails the stated ">2σ" test. At the four later times the smallest
adjacent difference is 8.6, 11.4, 13.3 and 10.8$\sigma$, so "8σ or more
thereafter" stands.

The end-to-end figure was added because it is both correct and much stronger:
across the five bins the difference is 9.0, 57, 62, 87 and 39$\sigma$. The
error bars behind this are plotted in Figure 3c, so a referee could have
checked the original claim.

### 1.2 "completely" vs "essentially" at 31.9 °C

Main text and SI disagreed about the same measurement.

- **Was (main):** "at 31.9 °C the correlation decays **completely** to its baseline"
- **Now (main):** "at 31.9 °C the correlation decays **essentially** to its baseline"

The SI has always said "decays essentially to the baseline". The weaker word is
the defensible one and is now used in both places.

### 1.3 Preposition

- **Was:** "the fraction of dynamics arrested **on** the experimental time window"
- **Now:** "the fraction of dynamics arrested **within** the experimental time window"

Everywhere else the paper writes "arrested within the accessible delay range"
or "arrested on that time scale".

---

## 2. Strengthening the reversibility claim

The title leads with reversible β-sheet-like ordering, but the seven-cycle test
of Figure S6 was run at 8-ID-I, which has no WAXS detector and reaches only
$Q \approx 0.033$ Å$^{-1}$. It therefore establishes reversibility of the
mesoscale structure, not of the 0.72 and 1.38 Å$^{-1}$ peaks. The evidence for
the peaks is at 12-ID-B, and it is stronger than the draft claimed: the
Measurement aliquot, on its own, had already been through more than 10 cycles.

**Added** after the peak positions:

> The Measurement aliquot alone makes the point for the high-$Q$ features:
> after more than 10 cycles it still develops both peaks at 30 °C and shows
> neither at 10 °C.

This converts a between-aliquot comparison into a within-aliquot demonstration
without adding any new measurement.

---

## 3. Symbol consistency (ACS Appendix 2)

The guide requires that symbols in the artwork match those in the text.

| | Was | Now |
|---|---|---|
| Scattering vector, text | $q$ | $Q$ (176 occurrences in main, SI, cover letter) |
| Scattering vector, figures | $Q$ | unchanged |
| Delay axis, Figures 3b, S4a, S5b, S6c, S10 | $\tau$ | $\Delta t$ |
| Delay symbol, Eq. 1 and Figure 1 | $\Delta t$ | unchanged |

$\tau$ was doing double duty: the text uses it for the fitted relaxation times
$\tau_\mathrm{fast}$ and $\tau_\mathrm{slow}$, so labelling the delay axis
$\tau$ made Figure 3b ambiguous against Figure 3c and Figure S9. Figure 1
already used $\Delta t$, so the two figures also disagreed with each other.

The conversion was scripted and restricted to math spans, so `\sqrt`,
`\qquad` and prose were untouched; no lowercase $q$ remains in any math span.

---

## 4. Wording

| Was | Now | Why |
|---|---|---|
| "a dispersed state at 10 °C" (abstract) | "a **molecularly dispersed** state at 10 °C" | Requested. Distinguishes the low-temperature state from the protein-poor phase of the assembled state, without claiming the narrow size distribution that "monodisperse" would assert and that a single Guinier fit at ~100 mg mL⁻¹ cannot support. |
| "…falls progressively with elapsed time, **so the central** result is a reversible structural transition…" | "…falls progressively with elapsed time. **The central result is therefore** a reversible structural transition…" | The "so" asserted that the XPCS decay implies the central result. It does not — the central result is the *combination* of the SAXS/WAXS reversibility and the XPCS arrest. Two sentences restore the original logic. |
| "coarsened, networklike morphology" (×2) | "coarsened, **network-like** morphology" | The summary paragraph already used the hyphenated form. |
| "0.25 and 0.34 Å⁻¹**,** (correlation lengths …)" | "0.25 and 0.34 Å⁻¹ (correlation lengths …)" | Stray comma before an opening parenthesis. |
| "peaks at $Q=1.38$ and $0.72$ Å⁻¹ compatible with β-sheet-like packing" | "peaks at $Q=1.38$ and $0.72$ Å⁻¹**,** compatible with…" | Without the comma the sentence reads as if the unit, not the peaks, is what is compatible. |

---

## 5. Supporting Information

**Added**, after the discussion of holding $\beta$ fixed:

> Nor does the choice drive the result: refitting every elapsed time with
> $\beta$ held at the extremes of the measured bin-to-bin spread, 0.1297 and
> 0.1328, shifts each fitted $f$ by less than 0.01 and leaves both
> monotonicities intact — $f$ decreasing with elapsed time at every $Q$, and
> increasing with $Q$ at every elapsed time.

$\beta$ is fixed and $f$ is the paper's only quantitative result, so a referee
will ask whether one is producing the other. The check was run; the numbers
above are its output.

**Section retitled.** "Public Availability of the Analysis Code and Data and
AI-Assisted Development" duplicated the preceding "Code and Reduced-Data
Availability" section. It is now "AI-Assisted Development". Both topics remain
listed in the SI contents and in the main-text SI description, which now ends
"…code and reduced-data availability; and AI-assisted development (PDF)."

---

## 6. Analysis code

`average_ranges.py` gained `spike_removal()`, a one-sided iterated modified
z-score (median ± 1.4826 × MAD, threshold 3) on the per-acquisition mean of
$I(Q)$ over 0.004–0.008 Å⁻¹ — the same band Figure S6b already reports its
repeatability in, so no new choice of range is introduced.

This **replaces** the eleven hand-picked `MANUAL_EXCLUDE` entries for B0083,
which are now selected automatically (ten of the same eleven; the cut also
takes frame 9 and keeps frame 47). `MANUAL_EXCLUDE` is now empty.

The cut removes 35 of the 1708 acquisitions surviving the cross-correlation
pass (2.0 %), applied to every group on identical terms:

- seven 6 °C groups: 4, 2, 1, 11, 4, 1, 1
- two buffer groups: 2, 7
- isothermal series: 2 — which move that group's absolute-scale coefficient by
  0.03 % and leave every fitted $g_2$ parameter unchanged
- fourteen high-temperature windows: none

Effect on Figure S6: the 6 °C cycle-to-cycle RSD falls from 13 % to 8.6 %
(quoted as 9 %). No fitting model, fit range, weighting, uncertainty
propagation, calibration or background subtraction was changed anywhere.

Also: `nexus_read.py` added (shared HDF readers, so all four 8-ID figure
scripts provably read the files the same way — it is imported by `saxpcs.py`
and was previously untracked, so the repository could not run without it);
`q_log_ticks()` in `acs_style.py`; panel-label offset raised to clear the
topmost y tick label.

---

## 7. Numbers re-verified, not changed

Re-derived by re-running the scripts against the committed data:

| Quantity | Script output | In the paper |
|---|---|---|
| 12-ID-B low-$Q$ slopes, 0.012–0.040 Å⁻¹ | −3.167, −3.146 | −3.17, −3.15 |
| 12-ID-B slopes, 0.012–0.033 Å⁻¹ (SI) | −3.091, −3.068 | −3.09, −3.07 |
| WAXS peak 1 | 0.7169 Å⁻¹, FWHM 0.337 | 0.72 Å⁻¹, 0.34 |
| WAXS peak 2 | 1.3829 Å⁻¹, FWHM 0.255 | 1.38 Å⁻¹, 0.25 |
| Guinier | $R_g = 22.649 \pm 0.502$ Å, $Q_\mathrm{max}R_g = 1.2457$ | 22.65 ± 0.50, 1.246 |
| $f$ at $t_w=5039$ s | 0.415–0.711 | 0.42–0.71 |
| $f$ at $t_w=7863$ s | 0.000–0.127 | 0.00–0.13 |
| $\gamma_\mathrm{fast}$ | −1.994 to −2.419 | −2.0 to −2.4 |
| reduced $\chi^2$ | 0.78–1.89 | 0.78–1.89 |
| $p_\mathrm{fast}$, $p_\mathrm{slow}$ | 0.355–0.703, 0.255–0.565 | 0.36–0.70, 0.26–0.57 |
| 6 °C RSD | 8.6 % | 9 % |
| 31.9 °C RSD, and detrended | 32.6 %, 7.1 % | 33 %, 7 % |
| 33.8 °C RSD, and detrended | 4.8 %, 1.4 % | 5 %, 1 % |
| Group sizes / durations | 63–200 acquisitions, 6.6–15.5 min | 63–200, 7–16 min |
| Frame time / max delay (from the HDF files) | 20 µs, 50 kHz, 1.6384 s | 20 µs, 50 kHz, 1.6 s |

Figure generation is deterministic: two consecutive runs of a figure script
produce PDFs that differ only in `/CreationDate`.

---

## 8. Kept as written

- **"The temperature-ramp data of Figure S6c show…"** — kept plural. "Data" is
  plural in formal scientific usage and ACS house style follows that.
- **Cover letter says "a dispersed state", not "molecularly dispersed"** — the
  extra word pushed the signature line onto a page of its own. The letter is
  back to two pages. The abstract keeps the new wording.
- **Figure 1's $g_2$ schematic** — the red curve is drawn parallel to the
  tilted $\Delta t$ axis, so the apparent rise is the panel's projection
  geometry, not a rising $g_2$. No change.
- **63–200 acquisitions per group** — correct: the range spans all six reported
  times including the 200-acquisition group at $t_w = 0$ shown in Figure 3a.
  The five analyzed XPCS groups are 63–150.
- **12-ID-B beamline citations** are papers that used the beamline rather than
  an instrument paper. There is no instrument paper; this is normal practice.
