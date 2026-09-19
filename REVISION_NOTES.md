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

---

## 9. Response to the 25 PDF comments

All 25 comments bound one-to-one to the 25 highlighted spans, in document
order. Line numbers are those of the commented draft.

### Page 2

**C1 · L15 · "Double-check all references in the Introduction section."**
Done for the whole bibliography, not only the Introduction, against the live
CrossRef record:

```
python3 tools/check_refs.py
```

**41 cited references; 40 journal articles all verified, 0 discrepancies.** The
41st is the 12-ID-B beamline page, a `@misc` web citation with no CrossRef
record. Each entry is looked up by DOI where it has one and by title otherwise,
and its year, volume, first page, journal and title are compared with the
publisher record.

Two genuine errors found and fixed — both titles that spelled out a Greek
letter the published title sets as a symbol:

| Key | Was | Now |
|---|---|---|
| `Serpell2000-pnas` | "…Synthetic **Alpha**-Synuclein Filaments Shows Amyloid-Like Cross-**Beta** Conformation" | "…Synthetic **$\alpha$**-Synuclein Filaments Shows Amyloid-Like Cross-**$\beta$** Conformation" |
| `Sawaya2007-nature` | "Atomic Structures of Amyloid Cross-**Beta** Spines…" | "Atomic Structures of Amyloid Cross-**$\beta$** Spines…" |

Both now render as α and β in the reference list, matching the published titles
and the β used throughout the paper. `Pan1993-pnas` keeps "Alpha-Helices" and
"Beta-Sheets" because the 1993 PNAS title really does spell them out.

Also added: a DOI for `Tanaka2000-qq` (10.1088/0953-8984/12/15/201). Nothing
was wrong with the entry — J. Phys.: Condens. Matter **12**, R207–R264 (2000)
is correct — but "Viscoelastic phase separation" is too generic a title for a
search to rank reliably, and the DOI pins the lookup for future runs.

Three journal names the first pass flagged are *not* errors and the checker no
longer reports them: "Annu. Rev. Biomed. Eng." and "Methods Mol. Biol." (×2)
are correct ACS abbreviations.

### Page 3

**C2 · L54 · "The observation of the paper is important. Break it into a new sentence."**
- **Was:** "…falls progressively with elapsed time, so the central result is a reversible structural transition coexisting with…"
- **Now:** "…falls progressively with elapsed time. The central result is therefore a reversible structural transition coexisting with…"

### Page 4

**C3 · L77 · "That's not the same with what the figure shows. Should we show the fitting in Fig. 2, or at least draw Q^-3.2 as a guide to the eye?"**
The fit *was* already drawn — a black dashed line over 0.012–0.040 Å⁻¹ — but it
runs straight through the 30 °C data it was fitted to and was invisible against
the markers at print size. It is now **red dashed with `zorder=5`**, matching
the colour this figure already uses for a fit (the two Gaussian peak fits in the
inset). Position unchanged, so it still shows the quality of the fit rather than
being an offset guide. The residual mismatch is rounding only: the in-figure
label prints `slope` to one decimal (−3.2) while the text quotes −3.17 and
−3.15.

**C4 · L78 · "Break into new sentence. The high Q peaks are important."**
- **Was:** "…for the Measurement aliquot) **and the same pair of high-$Q$ features, at** $Q=1.38$…"
- **Now:** "…for the Measurement aliquot). **Both also develop the same pair of high-$Q$ features, at** $Q=1.38$…"

**C5 · L87 · "I need to read this part to better understand how the error falls back to 7%. Honestly I think it's Ok even if the error is not 7%…"**
**Kept as written.** The detrending is what explains the one number in that
sentence that looks bad — the 33 % spread at the mid-ramp state. Dropping the
explanation would leave a 33 % sitting unexplained in the main text next to 9 %
and 5 %, which invites exactly the referee question the sentence forecloses.
The mechanism is in SI Section 6.2: the seven cycles cross acquisitions 241–250
between 31.62 and 32.48 °C, ln *I* tracks that window temperature at
*r* = 0.974 with d ln *I*/d*T* = 0.95 °C⁻¹, and removing the trend takes the
spread from 32.6 % to 7.1 %. Say the word and I will cut it.

### Page 5

**C6 · Figure 1 caption · "Yeah this definitely needs to be shorter … Focus more on what SAXS and SA-XPCS measures."**
Rewritten, **244 → 184 words**. Dropped: the explanation of the microfocused
beam, the Peltier/multi-zone thermal path (both already in the Experimental
Section), the beam-footprint bar, and the "one repeat unit … extended length"
gloss. Added, as the new centre of the caption:

> Averaging them gives the SAXS intensity $I(Q)$: how much structure the sample
> has at each length scale $2\pi/Q$, which grows at low $Q$ on heating.
> Correlating them in time gives $g_2(Q,\Delta t)$: how fast the structure at
> that same length scale rearranges, a decay meaning the density fluctuations at
> $2\pi/Q$ are still moving and a flat $g_2$ that they are frozen over the
> accessible delay.

The "illustrative, not data" and "conceptual illustration, not a real-space
observation" disclaimers are kept verbatim.

### Page 7

**C7 · L98 · "Need more explanation … is it really necessary over such a smoking gun evidence?"**
The hedge is kept — concurrence in two $Q$ ranges genuinely cannot identify the
contacts as the junctions — but it is no longer a bare assertion. **Added:**
"…: the same two signatures would appear if the ordered contacts sat inside the
protein-rich domains rather than at the links between them, and the two $Q$
ranges cannot distinguish those cases." A referee who sees *why* the caveat is
there reads it as rigour; one who sees only the caveat reads it as doubt.

**C8 · L104 · "The slow ramp protocol should be spelled out more clearly here."**
- **Was:** "…was brought to 30 °C on the slow ramp given in the Experimental Section."
- **Now:** "…was held at 6 °C for more than 30 min, warmed to 20 °C at 5 °C min⁻¹, and then brought to 30 °C at 1 °C min⁻¹, the slow final ramp limiting overshoot past the target."

**C9 · L105 · "Again, spell out where the RTD is."**
- **Now:** "the calibrated RTD, **embedded in the Peltier stage adjacent to the cell**, read 30 °C".

**C10 · L106 · "I won't call it 'first isothermal acquisition' … Change it to something like 'high temperature SAXS measurement shown in Fig. 2'."**
Phrase removed as asked, but **not** replaced with the suggested wording — see
§10 below. The sentence is about the 8-ID-I isothermal series, so it now reads
"…that moment is also the start of the first acquisition plotted in Figure 3."

**C11 · L107 · "I believe only one group is 63 measurements."**
Correct — only the last. **Was:** "a group of 63–200 consecutive 2 s
acquisitions". **Now:** "The groups are not all the same size. The five analyzed
for $g_2$ hold 150, 100, 100, 100 and 63 acquisitions … the last is short only
because the measurement ended there. The $t_w=0$ SAXS profile of Figure 3a
averages the first 200 acquisitions."

**C12 · L109 · "spell out the duration of the groups clearly."**
- **Was:** "a group spans 7–16 min" — which was also wrong, because it excluded
  the $t_w=0$ group.
- **Now:** "…and run 15.5, 10.4, 10.4, 10.7 and 6.7 min … [the $t_w=0$ group]
  over 20.8 min." Durations are from the recorded acquisition start times
  (cadence 6.2–6.5 s per 2 s acquisition, the rest being stage motion and
  readout).

**C13 · L116 · "Shouldn't this be 'shows'?"**
"Data … show" is correct — *data* is plural in ACS usage — so the original was
not an error. But the question is worth not provoking in a referee, so the
sentence was recast to make the subject singular and shorter:
- **Was:** "The temperature-ramp data of Figure S6c **show** the still-mobile limit directly:"
- **Now:** "Figure S6c **shows** the still-mobile limit directly:"

### Page 8

**C14 · L117 · "The legend in Fig. S6c shows 32 C. Please fix that."**
Real bug. `thermal_cycle.py` formatted the key rows with `:.0f`, rounding 31.9
and 33.8 to 32 and 34 while the text and caption said 31.9 and 33.8. Now `:.1f`
for the two ramp windows (6 °C stays `:.0f`). Figure regenerated; the key reads
6 °C / 31.9 °C / 33.8 °C.

**C15 · L129 · "1. … 'on the experimental time window' … makes the sentence very hard to read. 2. … should be broken into two sentences."**
- **Was:** "The fraction of fast dynamics therefore decreases while the fraction of dynamics arrested on the experimental time window increases, rather than all relaxation times shifting uniformly."
- **Now:** "The fraction of fast dynamics therefore decreases while the arrested fraction increases. What changes is the weight carried by each population; the relaxation times themselves do not simply shift together."

**C16 · L137 · "It's recommended to show in Eq. 2 that tau_fast and tau_slow are q dependent."**
Done — and `f` was given the same treatment, because leaving it as bare `f`
beside `τ(Q)` would imply the fast fraction is $Q$-independent, which is the one
thing the paper's main result contradicts. Eq. 2 and SI Eq. S4 now read
$f(Q)$, $\tau_\mathrm{fast}(Q)$, $\tau_\mathrm{slow}(Q)$;
$p_\mathrm{fast}$ and $p_\mathrm{slow}$ deliberately carry no argument, which
now encodes the global fit in the equation itself.

**C17 · L140 · "Are you referring to the global fitting of p_fast and p_slow? … emphasize that fewer fitting parameters reduce fit uncertainty."**
- **Now:** "The stretching exponents … carry no $Q$ argument because they are the
  shared parameters of a global fit: at each elapsed time the five measured $Q$
  bins were fitted simultaneously with both exponents common to all five.
  Sharing them reduces the free parameters from 25 to 17 and tightens the ones
  that remain." (25 = 5 bins × 5 parameters; 17 = 2 shared + 5 × 3.)

### Page 9

The whole page-9 paragraph was rebuilt around C18–C24 and split in two at
C22's request.

**C18 · L144 · "Spell out which value corresponds to which Q."**
- **Now:** "At $t_w=5039$ s it rises across the five bins from $f=0.42$ at
  $Q=0.00376$ Å⁻¹ (167 nm) to 0.71 at $Q=0.00827$ Å⁻¹ (76 nm); by $t_w=7863$ s
  the same bins give 0.00 to 0.13 … so $1-f$ has gone from 0.58 and 0.29 at
  those two $Q$ values to 1.00 and 0.87."

**C19 · L148 · "'an increasing fraction' over what? Q or t_w?"**
- **Now:** "**As $t_w$ grows**, an increasing fraction of the scattering ensemble
  is therefore arrested within the accessible delay range, and this holds at
  every length scale probed."

**C20 · L154 · "What is this statement doing here?"**
The clause "independent of elapsed time" was carrying an argument it never
stated. **Now:** "That ordering holds separately within each acquisition, so it
reflects length scale rather than the overall trend with $t_w$, and it ties the
growing low-$Q$ SAXS structure to the arrest fraction measured on the very same
frames. It is what coarsening predicts: the largest structural features become
the least mobile."

**C21 · L155 · "Add something like 'Note that' in the front."**
- **Now:** "**Note that** XPCS does not locate those dynamically arrested regions…"

**C22 · L157 · "It might be better to start a new paragraph … Right now everything is intertwined."**
New paragraph starts at "The auxiliary fit parameters behave consistently across
the series." The arrest result and the supporting fit parameters are now
separate paragraphs. A bulleted list was not used: Nano Letters Letters run as
continuous prose with no subsection structure.

**C23 · L159 · "Describe the trend briefly. How does it vary with t_w and Q?"**
- **Now:** "$\tau_\mathrm{fast}(Q)$ decreases with $Q$ at every analyzed time
  with an approximately diffusive-like exponent, $\gamma_\mathrm{fast}=-2.0$ to
  $-2.4$, and that exponent shows no systematic drift with $t_w$: the fast mode
  keeps the same $Q$ dependence while steadily losing amplitude."
  (Values in $t_w$ order: −2.29, −2.42, −2.20, −1.99, −2.35 — scattered, not
  monotonic.)

**C24 · L163 · "Highlight APS-U which provides the higher flux."**
- **Now:** "an independent control at **23–41 times the flux used here, a
  comparison the beam delivered by the APS Upgrade made possible**, shows no
  systematic flux-dependent change…"

### Page 11

**C25 · L172 · "Don't you want to add that the higher Q lost mobility faster …?"**
Added, but **with the direction reversed** — see §10.
- **Now:** "…the low-$Q$ assembly progressively loses mobility, **doing so
  fastest at the lowest $Q$ — the largest structures stop moving first, as a
  coarsening network should.**"

---

## 10. Two comments not implemented as written

**C10 — the suggested replacement points at the wrong dataset.** The sentence is
about the 8-ID-I isothermal series: it defines $t_w=0$ as the moment the RTD
reached 30 °C. Redirecting it to "high temperature SAXS measurement shown in
Fig. 2" would point the reader at the 12-ID-B capillary data, which is a
different sample, a different beamline and has no waiting-time axis at all —
the SI says so explicitly. The objectionable phrase is gone; the pointer now
goes to Figure 3, which is the isothermal series being defined.

**C25 — the direction is the other way round.** $f$ is the **fast** fraction and
it *increases* with $Q$, so the arrested fraction $1-f$ is largest at the
**lowest** $Q$. At $t_w=7863$ s, $1-f = 1.00$ at $Q=0.00376$ Å⁻¹ (167 nm) and
0.87 at $Q=0.00827$ Å⁻¹ (76 nm). Mobility is therefore lost fastest at the
largest length scale, not the smallest. Written as "the higher $Q$ lost mobility
faster" the sentence would contradict Figure 3c and the paragraph on page 9.
The coarsening reading you want is right and is now stated in the summary — it
just runs low-$Q$-first.

---

## 11. Updated counts

Abstract 129 words (limit 150). Main text **2782 words excluding captions**
(limit 3000), 3280 including them. Captions: Figure 1 184, Figure 2 119,
Figure 3 195. Main PDF is 19 pages; `make check` is clean.

---

## 12. A less painful way to send comments next time

Screenshots of the Comments pane arrive detached from the highlights, so the
notes and the marked text have to be matched up by hand. Three options, best
first.

**1. Export the comments and commit them (recommended).** In Acrobat: Comments
pane → the "…" menu → *Export All to Data File* → save as `.xfdf` under
`review/` in this repo, alongside the PDF you annotated. Then:

```
python3 tools/pdf_comments.py review/comments.xfdf manuscript/build/main.pdf
```

prints every note next to the exact words it marks, in document order:

```
--- C7  page 7
    marked : but does not by itself prove that the contacts are the network junctions
    comment: Need more explanation. I understand you are trying to sound safe...
```

`.xfdf` is XML, so the note text is directly readable, and the word boxes come
from `pdftotext -bbox-layout`. Standard library and poppler only — no PDF
package needed, none is installed on the beamline machines. The script is
tested on single-line and multi-line spans.

**2. Commit the annotated PDF itself.** Dropping `main_QZ_commented.pdf` into
`review/` also works and needs no export step; the annotations can be read off
it directly. Slightly less robust than the `.xfdf`, because Acrobat may store
the annotation dictionaries inside compressed object streams.

**3. Line numbers, no tooling at all.** `main.pdf` already carries line numbers
down the margin — the `lineno` package is loaded for exactly this. A plain file
like `review/comments.md` with

```
L54  break into a new sentence, the observation is the important part
L107 only one group is 63 -- spell the sizes out
L117 Fig. S6c key says 32 C
```

is unambiguous and takes no parsing. Good for a short pass; the `.xfdf` route
is better when there are twenty-five of them.

What does **not** work is pasting comment text without the anchors, or
screenshots of the Comments pane: both drop the link between each note and its
highlight, which is the part that costs the time.
