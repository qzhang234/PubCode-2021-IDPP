# Citation audit

Every reference cited in the manuscript or the Supporting Information: its DOI,
where it is cited, and what it is doing there.

41 references, 65 citation instances. All 40 journal articles were checked
against the live CrossRef record with `python3 tools/check_refs.py` — year,
volume, first page, journal and title all agree. The 41st is a beamline web
page with no CrossRef record. DOIs for the 21 entries whose BibTeX record
lacked one were resolved from CrossRef and are given below; they are not all in
`reference.bib`, which is fine for Fast Format but worth adding if the journal
asks for DOIs at revision.

Grouped by the job each reference does in the argument, following the order of
the manuscript.

---

## A. Amyloid as a pathological outcome

One reference per disease, supporting the opening list. None is doing more work
than naming the disease and the protein involved.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Masters *et al.*, *PNAS* **82**, 4245–4249 (1985) | [10.1073/pnas.82.12.4245](https://doi.org/10.1073/pnas.82.12.4245) | Intro, ¶1 | Identifies the Aβ plaque core protein — the Alzheimer's entry in the disease list. |
| Serpell *et al.*, *PNAS* **97**, 4897–4902 (2000) | [10.1073/pnas.97.9.4897](https://doi.org/10.1073/pnas.97.9.4897) | Intro, ¶1 | Fiber diffraction showing α-synuclein filaments are cross-β — the Parkinson's entry. Also the closest methodological precedent for reading β-sheet packing off a diffraction pattern, which is what Figure 2's inset does. |
| Cooper *et al.*, *PNAS* **84**, 8628–8632 (1987) | [10.1073/pnas.84.23.8628](https://doi.org/10.1073/pnas.84.23.8628) | Intro, ¶1 | Isolation of amylin from type II diabetic pancreas — the diabetes entry. |
| Pan *et al.*, *PNAS* **90**, 10962–10966 (1993) | [10.1073/pnas.90.23.10962](https://doi.org/10.1073/pnas.90.23.10962) | Intro, ¶1 | α-helix → β-sheet conversion in scrapie prion protein — the prion entry. |
| Radamaker *et al.*, *Nat. Commun.* **10**, 1103 (2019) | [10.1038/s41467-019-09032-0](https://doi.org/10.1038/s41467-019-09032-0) | Intro, ¶1 | Cryo-EM structure of an AL amyloid fibril — the light-chain amyloidosis entry, and the most recent structural example in the list. |

## B. Why amyloid-like order is not inherently pathological

These carry the pivot the introduction turns on: the same fold appears in
functional biology and in engineered materials.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Sawaya *et al.*, *Nature* **447**, 453–457 (2007) | [10.1038/nature05695](https://doi.org/10.1038/nature05695) | Intro ¶1; Results ¶4 | **Cited twice.** In the introduction it supports *why* many fibrils are kinetically persistent — the dry, tightly packed steric-zipper interface. In the results it is one of the two sources for reading 4.5 and 8.8 Å as strand-to-strand and sheet-to-sheet spacings. The contrast between a steric zipper and what we see is the point: our peaks are broad, so the ordering runs over a few repeats, not a crystallite. |
| Maji *et al.*, *Science* **325**, 328–332 (2009) | [10.1126/science.1173155](https://doi.org/10.1126/science.1173155) | Intro, ¶1 | Peptide hormones stored as functional amyloid in pituitary granules, which *disassemble* to release monomer. This is the load-bearing citation for "what differs is how long it lasts" — a natural precedent for reversible amyloid-like packing, which is the paper's thesis. |
| Chapman *et al.*, *Science* **295**, 851–855 (2002) | [10.1126/science.1067484](https://doi.org/10.1126/science.1067484) | Intro, ¶1 | *E. coli* curli as functional amyloid in biofilm matrix — the second functional-amyloid example, and a structural rather than storage role. |
| Knowles & Mezzenga, *Adv. Mater.* **28**, 6546–6561 (2016) | [10.1002/adma.201505961](https://doi.org/10.1002/adma.201505961) | Intro, ¶1 | Review of amyloid fibrils as building blocks for functional materials — establishes that the materials-design framing of the paper is an existing field, not a novelty claim. |

## C. Reversibility of amyloid-like contacts

The specific literature that makes "reversible β-sheet-like junction" a
credible hypothesis rather than a contradiction in terms.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Pálmadóttir *et al.*, *Biophys. Rev.* **6**, 011303 (2025) | [10.1063/5.0236947](https://doi.org/10.1063/5.0236947) | Intro ¶1 and ¶2 | **Cited twice.** The review that licenses the central framing: amyloid formation spans a range of reversibility rather than being a one-way commitment. The newest reference in the paper and the one a referee is most likely to check. |
| Kato *et al.*, *Cell* **149**, 753–767 (2012) | [10.1016/j.cell.2012.04.017](https://doi.org/10.1016/j.cell.2012.04.017) | Intro ¶2; Results ¶4 | **Cited twice.** Low-complexity domains forming dynamic, reversible cross-β fibers inside hydrogels. Directly analogous to the claim we make: β-sheet-like contacts that act as physical junctions and come apart again. |
| Hughes *et al.*, *Science* **359**, 698–701 (2018) | [10.1126/science.aan6398](https://doi.org/10.1126/science.aan6398) | Intro ¶2; Results ¶4 | **Cited twice.** Atomic structures of kinked β-sheets that assemble networks — the structural reason a hydrated, weakly packed interface can be labile where a steric zipper is not. Pairs against Sawaya as the counter-example. |
| Luo *et al.*, *Nat. Struct. Mol. Biol.* **25**, 341–346 (2018) | [10.1038/s41594-018-0050-8](https://doi.org/10.1038/s41594-018-0050-8) | Intro ¶2; Results ¶4 | **Cited twice.** FUS LC segments giving explicitly *reversible* amyloid fibrils, with the structural basis resolved. The closest published precedent for what we infer from the disappearance of the WAXS peaks on cooling. |

## D. ELP sequence design and LCST behaviour

Establishes that the platform is programmable, which is what makes the result a
materials-design statement rather than a one-off observation.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Meyer & Chilkoti, *Biomacromolecules* **5**, 846–851 (2004) | [10.1021/bm034215n](https://doi.org/10.1021/bm034215n) | Intro, ¶2 | Quantifies how chain length and concentration set the ELP transition temperature — the basis for "LCST behavior can be programmed". |
| Amiram *et al.*, *Nat. Mater.* **10**, 141–148 (2011) | [10.1038/nmat2942](https://doi.org/10.1038/nmat2942) | Intro ¶2; SI §1.2 | **Cited twice, for two different things.** In the introduction, recombinant control of ELP sequence and molecular weight. In the SI it is the actual method reference for OERCA, the technique used to build the (VPAVG)₃₀ gene. |
| Quiroz & Chilkoti, *Nat. Mater.* **14**, 1164–1171 (2015) | [10.1038/nmat4418](https://doi.org/10.1038/nmat4418) | Intro, ¶2 | Sequence heuristics for encoding phase behaviour in disordered protein polymers — supports the claim that phase behaviour is sequence-encoded and therefore transferable. |
| Varanko *et al.*, *Annu. Rev. Biomed. Eng.* **22**, 343–369 (2020) | [10.1146/annurev-bioeng-092419-061127](https://doi.org/10.1146/annurev-bioeng-092419-061127) | Intro ¶2 (×2) | Review, cited both for recombinant programmability and for the coacervate/assembly behaviour above the LCST. |
| Saha *et al.*, *Adv. Ther.* **3**, 1900164 (2020) | [10.1002/adtp.201900164](https://doi.org/10.1002/adtp.201900164) | Intro ¶2 (×2) | Review of ELP architectures from unimers to hierarchical self-assembly — supports "hierarchical assemblies depending on sequence, concentration, and thermal history", which is the regime Figure 2 shows. |

## E. ELP gels formed by arrested phase separation

The immediate prior art. The paper's novelty claim is measured against these,
so they are the references a referee will weigh hardest.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Glassman & Olsen, *Biomacromolecules* **16**, 3762–3773 (2015) | [10.1021/acs.biomac.5b01026](https://doi.org/10.1021/acs.biomac.5b01026) | Intro ¶2; Fig. 1 text; Fig. 1 caption; Results ¶4; Discussion ¶1 | **Cited five times — the most-used reference in the paper.** The closest prior system: alanine-containing ELPs forming stiff thermoresponsive gels by arrested phase separation. It supports the stated gap ("Building on that work, the complementary question…"), it is one of the two justifications for the Figure 1 cartoon, and it supplies the bicontinuous/fractal network morphology we read the *q*⁻³·² regime against. |
| Glassman *et al.*, *Biomacromolecules* **17**, 415–426 (2016) | [10.1021/acs.biomac.5b01210](https://doi.org/10.1021/acs.biomac.5b01210) | Intro ¶2; Results ¶4; Discussion ¶1 | **Cited three times.** The follow-up: toughened arrested ELP networks. Supports the same three claims as above and extends them to mechanical behaviour. |
| Sing *et al.*, *Soft Matter* **13**, 8511–8524 (2017) | [10.1039/c7sm00638a](https://doi.org/10.1039/c7sm00638a) | Intro, ¶2 | Structure and rheology of associative protein hydrogels — supports that assembly above the LCST depends on sequence, concentration and thermal history. |
| Sing *et al.*, *Macromolecules* **51**, 2951–2960 (2018) | [10.1021/acs.macromol.8b00002](https://doi.org/10.1021/acs.macromol.8b00002) | Discussion, ¶1 | Long stress-relaxation times and gel-like elasticity in related ELPs — the rheological counterpart to the slow mode we resolve. This is the reference that makes our *τ*_slow physically plausible rather than a fitting artefact. |

## F. Aging and heterogeneous dynamics measured by XPCS

Establishes both that the measurement is standard practice and that the
two-step decay we see has precedent.

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Bahadur *et al.*, *J. Chem. Phys.* **151**, 104902 (2019) | [10.1063/1.5111521](https://doi.org/10.1063/1.5111521) | Intro ¶3; Results ¶7; SI §6.4 | **Cited three times.** XPCS plus rheology on thermo-reversible nanoparticle gels. Supports physical aging as a defined phenomenon, the two-step decay as a signature of heterogeneous arrest, and stretched exponents as evidence of a broad relaxation-time distribution. |
| Chen *et al.*, *J. Chem. Phys.* **158**, 024906 (2023) | [10.1063/5.0126432](https://doi.org/10.1063/5.0126432) | Intro, ¶3 | Memory in aging colloidal gels — supports that microscopic relaxation keeps evolving while time-averaged structure barely changes, the premise for needing a dynamical measurement at all. |
| Girelli *et al.*, *Phys. Rev. Lett.* **126**, 138004 (2021) | [10.1103/PhysRevLett.126.138004](https://doi.org/10.1103/PhysRevLett.126.138004) | Intro ¶3; Discussion ¶1 | **Cited twice.** XPCS of protein liquid–liquid phase separation and domain coarsening — the closest published analogue of this experiment, and part of the support for reading our result as network-forming arrested phase separation. |
| Begam *et al.*, *Phys. Rev. Lett.* **126**, 098001 (2021) | [10.1103/PhysRevLett.126.098001](https://doi.org/10.1103/PhysRevLett.126.098001) | Intro ¶3; Results ¶7; Discussion ¶1; SI §6.4 | **Cited four times.** Coherent X-ray study of egg-white gelation: network formation followed by heterogeneous dynamics. The closest precedent for the *specific* observable we report — a growing slow fraction during protein gelation. |
| Chushkin *et al.*, *Phys. Rev. Lett.* **129**, 238001 (2022) | [10.1103/PhysRevLett.129.238001](https://doi.org/10.1103/PhysRevLett.129.238001) | Intro, ¶3 | Cage relaxation in concentrated protein solutions by XPCS — establishes that XPCS resolves microscopic dynamics in crowded protein systems, which is the regime of our ~100 mg mL⁻¹ sample. |
| Chen *et al.*, *Phys. Rev. E* **102**, 042619 (2020) | [10.1103/PhysRevE.102.042619](https://doi.org/10.1103/PhysRevE.102.042619) | Results, ¶7 | Microscopic ergodicity breaking in glass-forming nanoclay suspensions — supports the two-step relaxation with a growing slow component as a general signature of kinetic arrest, in a non-protein system. |

## G. Theory and scattering interpretation

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Tanaka, *J. Phys.: Condens. Matter* **12**, R207–R264 (2000) | [10.1088/0953-8984/12/15/201](https://doi.org/10.1088/0953-8984/12/15/201) | Fig. 1 text; Fig. 1 caption; Discussion ¶1 | **Cited three times.** The viscoelastic phase-separation framework. It is the source of the dynamic-asymmetry picture drawn in Figure 1 stages i–iv, and the theoretical frame for reading a connected dense phase with suppressed relaxation. Flagged in both the text and the caption as motivating an illustration, not as something the data demonstrate. |
| Teixeira, *J. Appl. Crystallogr.* **21**, 781–785 (1988) | [10.1107/S0021889888000263](https://doi.org/10.1107/S0021889888000263) | Results, ¶4 | The standard reference for reading power-law small-angle scattering as self-similar density heterogeneity. Supports the deliberately careful phrasing that *q*⁻³·² means fractal-like heterogeneity over the probed interval, not one resolved structural scale. |
| Lu *et al.*, *Nature* **453**, 499–503 (2008) | [10.1038/nature06931](https://doi.org/10.1038/nature06931) | Discussion, ¶1 | Gelation of short-range-attractive particles as arrested spinodal decomposition — the colloidal-physics reference for the mechanism we invoke, generalising it beyond protein systems. |

## H. Cross-β spacings

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Sunde *et al.*, *J. Mol. Biol.* **273**, 729–739 (1997) | [10.1006/jmbi.1997.1348](https://doi.org/10.1006/jmbi.1997.1348) | Results, ¶4 | The canonical synchrotron X-ray diffraction assignment of the ~4.7 Å strand-to-strand and ~10 Å sheet-to-sheet spacings of the cross-β core. With Sawaya 2007 it is the entire basis for calling our 1.38 and 0.72 Å⁻¹ peaks β-sheet-like — which is why the limitations paragraph says the assignment rests on peak positions and not on spectroscopy. |

## I. Beamlines and instrumentation

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| APS, 12-ID-B beamline page (accessed 22 Jul 2026) | [aps.anl.gov/Sector-12/12-ID-B](https://www.aps.anl.gov/Sector-12/12-ID-B) | Results, ¶3 | Beamline and detector configuration for the SAXS/WAXS measurement. The only non-journal citation; there is no 12-ID-B instrument paper. |
| Xu *et al.*, *Nat. Commun.* **13**, 126 (2022) | [10.1038/s41467-021-27742-2](https://doi.org/10.1038/s41467-021-27742-2) | Results ¶3; SI §2.1 | Cited as a published example of the simultaneous SAXS/WAXS configuration at 12-ID-B, standing in for an instrument paper. |
| Zheng *et al.*, *ACS Nano* **16**, 4813–4822 (2022) | [10.1021/acsnano.2c00161](https://doi.org/10.1021/acsnano.2c00161) | Results ¶3; SI §2.1 | Second published use of the same 12-ID-B configuration, by beamline staff who are coauthors here. |
| Zhang *et al.*, *J. Synchrotron Radiat.* **28**, 259–265 (2021) | [10.1107/S1600577520014319](https://doi.org/10.1107/S1600577520014319) | SI §2.2 | The 20 µs-resolved XPCS capability on a 500k-pixel detector — the method reference for the 50 kHz acquisition and the multi-tau delay list reaching 1.638 s. |
| Sheyfer *et al.*, *Phys. Rev. Lett.* **125**, 125504 (2020) | [10.1103/PhysRevLett.125.125504](https://doi.org/10.1103/PhysRevLett.125.125504) | SI §3.2 | Source of the absolute-intensity conversion, in particular the per-pixel solid angle ΔΩ = (*d*_pix/*D*)², used to put Figure 3a on an absolute cross-section scale. |
| Ozgulbas *et al.*, *Light Sci. Appl.* **12**, 196 (2023) | [10.1038/s41377-023-01233-z](https://doi.org/10.1038/s41377-023-01233-z) | SI §2.2 | The multi-zone Peltier sample environment and its calibrated RTD readout — the provenance of the temperature control on which the whole *t*_w axis depends. |

## J. Molecular-biology methods (SI only)

| Reference | DOI | Cited at | Connection |
|---|---|---|---|
| Aslanidis & de Jong, *Nucleic Acids Res.* **18**, 6069–6074 (1990) | [10.1093/nar/18.20.6069](https://doi.org/10.1093/nar/18.20.6069) | SI §1.2 | Original ligation-independent cloning method, used to introduce the LIC adapters. |
| Eschenfeldt *et al.*, *Methods Mol. Biol.* **498**, 105–115 (2009) | [10.1007/978-1-59745-196-3_7](https://doi.org/10.1007/978-1-59745-196-3_7) | SI §1.2 | The pMCSG LIC vector family, of which pMCSG53 — the vector actually used — is a member. |
| Makowska-Grzyska *et al.*, *Methods Mol. Biol.* **1140**, 89–105 (2014) | [10.1007/978-1-4939-0354-2_7](https://doi.org/10.1007/978-1-4939-0354-2_7) | SI §1.3 | The *E. coli* expression and lysis protocol followed for production and purification. |

---

## Observations from the audit

- **Citation load is concentrated where it should be.** Glassman 2015 (5×),
  Begam 2021 (4×), Tanaka 2000, Glassman 2016 and Bahadur 2019 (3× each) are
  the most-cited — the immediate prior art in ELP gels, in XPCS of protein
  gelation, and the theory frame. That matches what the paper claims to build on.
- **Every claim in the four-limitation paragraph has its counter-citation
  nearby.** The β-sheet assignment rests on Sunde 1997 and Sawaya 2007 and is
  explicitly flagged as positional, not spectroscopic.
- **Nothing is cited for a claim it does not make.** The two places where the
  literature is used for something the data cannot show — the Figure 1 cartoon
  (Tanaka, Glassman 2015) and the junction identity (Kato, Hughes, Luo) — both
  carry explicit hedges in the text.
- **Single-use references are all doing real work**; there is no padding. The
  five disease citations in the opening list are the only place where one
  reference stands for one fact, and that is the conventional form for such a list.
- **One page field was wrong and is now fixed.** `Ozgulbas2023-oz` carried
  `pages = "1--10"`, a reference-manager placeholder; *Light: Science &
  Applications* numbers articles, and the real locator is **196**. It survived
  the first CrossRef pass because that record leaves `page` empty and puts the
  number in `article-number`, so the page comparison was skipped. Twelve
  references are numbered that way; the other eleven were all correct.
  `tools/check_refs.py` now falls back to `article-number`.
- **19 of the 40 DOIs are not in `reference.bib`.** Harmless for Fast Format,
  which accepts any complete reference style, but worth adding at revision —
  the resolved DOIs are in the tables above. (DOIs were added to
  `Tanaka2000-qq` and `Ozgulbas2023-oz` during this audit, to pin future
  lookups.)
