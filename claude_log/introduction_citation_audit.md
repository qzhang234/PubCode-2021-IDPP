# Introduction citation audit and figure review

**Manuscript:** *Reversible beta-Sheet-like Ordering and XPCS-Resolved Dynamic Arrest during Elastin-like Polypeptide Phase Separation* (Nano Letters submission)  

**Date:** 2026-08-23  

**Method:** 28 agents. Each Introduction reference was checked against Crossref, PubMed/Europe PMC, OpenAlex and publisher full text; every citation flagged as problematic was then re-checked by a second, adversarial agent instructed to refute the first. Three further agents audited figure/data explanation, cross-reference integrity, and whether the results support the argument. Findings that would change the paper's claims were independently re-verified by hand against the deposited data and code; those verifications are recorded in the last section.


---

## 1. Summary

All 26 distinct references cited in the Introduction exist and were located. No citation is fabricated and no reference is off-topic.


| # | Ref | Para | Claim supported? | Bib data | Action |
|---|---|---|---|---|---|
| 1 | `Chapman2002-science` | 1 | SUPPORTS | OK | none |
| 2 | `Cooper1987-pnas` | 1 | SUPPORTS | OK | none |
| 3 | `Glenner1971-science` | 1 | SUPPORTS | **errors** | author names fixed |
| 4 | `Maji2009-science` | 1 | SUPPORTS | **errors** | author names fixed |
| 5 | `Masters1985-pnas` | 1 | SUPPORTS | OK | none |
| 6 | `Palmadottir2025-bpr` | 1 | SUPPORTS | OK | none |
| 7 | `Pan1993-pnas` | 1 | SUPPORTS | OK | none |
| 8 | `Sawaya2007-nature` | 1 | SUPPORTS | OK | none |
| 9 | `Serpell2000-pnas` | 1 | SUPPORTS | OK | none |
| 10 | `Amiram2011-la` | 2 | SUPPORTS | OK | none |
| 11 | `Glassman2015-cw` | 2 | SUPPORTS | OK | none |
| 12 | `Glassman2016-pq` | 2 | SUPPORTS | OK | none |
| 13 | `Hughes2018-science` | 2 | SUPPORTS | **errors** | author name fixed |
| 14 | `Kato2012-cell` | 2 | PARTIAL | OK | kept (see notes) |
| 15 | `Luo2018-nsmb` | 2 | SUPPORTS | OK | none |
| 16 | `Meyer2004-zb` | 2 | SUPPORTS | OK | none |
| 17 | `Quiroz2015-mh` | 2 | SUPPORTS | OK | none |
| 18 | `Saha2020-jg` | 2 | SUPPORTS | OK | none |
| 19 | `Sing2017-zq` | 2 | SUPPORTS | OK | none |
| 20 | `Sing2018-er` | 2 | PARTIAL | OK | **removed from this claim** |
| 21 | `Varanko2020-ek` | 2 | SUPPORTS | OK | none |
| 22 | `Bahadur2019-ln` | 3 | SUPPORTS | OK | none |
| 23 | `Begam2021-fh` | 3 | SUPPORTS | OK | none |
| 24 | `Chen2023-ho` | 3 | SUPPORTS | OK | none |
| 25 | `Chushkin2022-bw` | 3 | SUPPORTS | OK | none |
| 26 | `Girelli2021-zz` | 3 | SUPPORTS | OK | none |

### Changes made

1. **`reference.bib`** - five wrong author given names corrected: Glenner 1971 (`Harada, Masahiro`->`Minoru`, `Isersky, Charles`->`Chaviva`), Maji 2009 (`Perrin, Michele H.`->`Marilyn H.`, `Jessberger, Sibylle`->`Sebastian`, `Nilsson, Karl P. R.`->`K. Peter R.`), Hughes 2018 (`Chong, Lin`->`Lisa`).  
   *Severity note:* achemso renders every author as initials only, and in each case the wrong given name shares its initial. **None of these errors ever reached the printed page** - this is database hygiene, not a submission defect.

2. **`main.tex`** - the clause "related assemblies are increasingly used as mechanically robust biomaterials" carried no citation of its own (Maji covers hormone storage, Chapman covers biofilm ECM; neither is about biomaterials). Added Knowles & Mezzenga, *Adv. Mater.* **2016**, 28, 6546-6561, doi:10.1002/adma.201505961 - verified against Crossref and added to the .bib as `Knowles2016-am`.

3. **`main.tex`** - `Sing2018-er` removed from "alanine-containing ELPs can form stiff, thermoresponsive gels through arrested phase separation". Sing 2018 studies *dual-associating ABA triblock* fusion proteins whose gelation comes from coiled-coil midblock association plus equilibrium endblock aggregation, and it *cites Glassman 2015* for arrested phase separation rather than demonstrating it. Glassman 2015 and 2016 each independently establish all three attributes of the claim (alanine XPAVG sequence, MPa-scale stiffness, arrested phase separation), so the sentence stands on them alone.  
   *Optional:* Sing 2018 is now uncited. It could be restored on a claim it does demonstrate, e.g. "Substituting alanine for glycine at the third position slows end-block association dynamics and alters deformation behaviour."


---

## 2. Per-citation detail: where in each cited paper the support lies

This is the record requested - for every Introduction reference, the specific section, figure or paragraph of the *cited* work that carries the claim.


### Introduction paragraph 1


#### `Chapman2002-science` - SUPPORTS

**Claim it is cited for:** ... and extracellular-matrix formation ...


**Supporting passage in the cited paper:**

> I independently scraped the complete full text from PMC2838482 (NIHMS author manuscript, publicly readable) rather than relying on the first pass. The assigned sub-clause -- functional amyloids mediate extracellular-matrix formation -- is supported by the conjunction of the paper's own central result and its own first-hand imaging. (1) THE AMYLOID IDENTITY IS THIS PAPER'S OWN RESULT. Abstract, verbatim from PubMed: "Biochemical, biophysical, and imaging analyses revealed that fibers produced by Escherichia coli called curli were amyloid." Results, Fig. 2, verbatim: "Circular dichroism (CD) analysis indicated that these fibers were rich in beta-sheet secondary structure with a minimum peak at ~218 nm (Fig. 2A)" and "Like other amyloid fibers, S6 curli induced a spectral change of a 10 uM Congo red (CR) solution..." (2) THE MATRIX MORPHOLOGY IS ALSO THIS PAPER'S OWN OBSERVATION, not inherited background. Results, describing the authors' own deep-etch EM in Fig. 1A, verbatim from my scrape: "Under high-resolution EM, curli appeared as a tangled and amorphous matrix surrounding the bacteria (Fig. 1A)." (3) THE EXTRACELLULAR, MATRIX-FORMING FUNCTION is stated in the paper's opening paragraph, verbatim: "Curli are a class of highly aggregated, extracellular fibers expressed by Escherichia and Salmonella spp. that are involved in the colonization of inert surfaces and biofilm formation (1, 2) and mediate binding to a variety of host proteins (3-5)." (4) Concluding paragraph, verbatim: "Our demonstration that E. coli can produce extracellular amyloid-like fibers increases the recognized functional repertoire of amyloid fibers and provides a useful model system to study their formation." WHAT I CONFIRM IS ABSENT: the exact adjacent phrase "extracellular matrix" does not occur (the two words appear in the paper describing curli, but in different sentences); and a keyword scan returns 0 hits for "biomaterial", "mechanic", "stiff", "elastic", "robust", "tough", "engineer". The only physical-property remark is qualitative: the CsgA-his solutions "became opaque and noticeably viscous."


#### `Cooper1987-pnas` - SUPPORTS

**Claim it is cited for:** ... type II diabetes


**Supporting passage in the cited paper:**

> Support is concentrated in three places, all read directly from the full text (Europe PMC scanned PDF of PMC299599).
> 
> (1) INTRODUCTION, opening two sentences, p. 8628 left column - the disease association plus the fibrillar nature of the deposit:
> "The occurrence of amyloid in the islets of Langerhans is a major pathological feature of type 2 diabetes."
> and, three sentences later:
> "This was later confirmed by alkaline Congo red staining (4) and by the finding of fibrillar structure on electron microscopy (5)."
> 
> (2) MATERIALS AND METHODS, subsection "Human Subjects," p. 8628 right column - the amyloid is operationally confirmed by the classic Congo red criterion and is found only in the diabetic pancreases:
> "Islet amyloid was detected in tissue fixed in 150 mM NaCl/10% formalin by light microscopy after hematoxylin and eosin staining and was confirmed after staining with alkaline Congo red by the demonstration of green birefringence by microscopy under polarized light. Only the pancreases from diabetic patients contained amyloid."
> 
> (3) INTRODUCTION final paragraph and DISCUSSION first paragraph, pp. 8628 and 8630 - the specific precursor protein (DAP, now known as amylin / islet amyloid polypeptide, IAPP) and its diabetes-specific occurrence:
> "Quantitative studies of protein yields indicated DAP is the major protein constituent of islet amyloid."
> and
> "neither particulate amyloid nor DAP was detectable in any of six amyloid-negative pancreases from nondiabetic subjects. In another study, we have shown that pancreatic islet amyloid deposits were found in 22 of 24 type 2 diabetic subjects and were not present in 10 age-matched control subjects (25)."
> 
> Also directly relevant to the manuscript's beta-sheet framing, RESULTS, final sentence before DISCUSSION, p. 8630 left column:
> "These studies suggest a tendency to [beta]-sheet formation in the middle of the molecule where DAP is hydropathic."
> (The source is a 1987 scan; OCR renders the Greek beta as "8" here and as "13" in the abstract's "beta-pleated sheets." I have substituted "[beta]" and flag it rather than silently correcting. All other quotes above are verbatim from the OCR text and contain no OCR ambiguity.)


#### `Glenner1971-science` - SUPPORTS

**Claim it is cited for:** ... systemic light-chain amyloidosis


**Supporting passage in the cited paper:**

> Abstract, which for this two-page Science report is also the whole thesis. Verbatim, retrieved independently from BOTH Crossref's JATS abstract deposit and Europe PMC (the Science full text is genuinely closed: www.science.org returns HTTP 403, Unpaywall reports is_oa=false / oa_status=closed with zero OA locations, no PMC deposit):
> 
> "The sequences of the 35 and 36 amino-terminal amino acids of two purified amyloid fibril proteins have been determined. Results indicate that these two proteins are derived from homogeneous immunoglobulin light chains of variable region subgroup V(kappaI). The relation between amyloidosis and immunoglobulins is thus more firmly established."
> 
> The operative clauses are "these two proteins are derived from homogeneous immunoglobulin light chains" and "The relation between amyloidosis and immunoglobulins is thus more firmly established." The title is itself the claim: "Proof of Homology with Immunoglobulin Light Chains by Sequence Analyses."
> 
> The underlying evidence is the residue-by-residue N-terminal alignment of the two purified fibril proteins against V(kappa)I light-chain prototypes, in the body of the report. I could not read the body, so I cannot name a figure or table number and have no verbatim body quote.
> 
> Corroborating indexing (MEDLINE, PMID 4102463): major MeSH headings include "Amyloidosis / etiology / immunology" alongside "Amyloid / analysis" and "gamma-Globulins / analysis" - NLM indexed this paper as establishing the etiology of amyloidosis.
> 
> On the "systemic" qualifier specifically, I obtained a source the first pass did not have in full text. Kisilevsky, Raimondi & Bellotti, Front. Mol. Biosci. 2016, 3, 17 (PMC4860540), full text retrieved: "In 1970 the first amyloid protein was isolated from natural deposits occurring in a patient affected by multiple myeloma and amyloidosis. This amyloid protein was composed mainly of a fragment of a monoclonal light chain whose primary structure was identical to the variable region of the monoclonal light chains isolated from the patient's urine (Glenner et al., 1970, 1971a, c)." The same review labels this disease category "primary (now AL)". Myeloma-associated / plasma-cell-dyscrasia amyloidosis is systemic by definition, so the source material behind Glenner 1971 is systemic light-chain amyloidosis, not localized deposits.


**Bibliographic corrections:** UPHELD (I tried to refute this and could not; the evidence against the .bib is now stronger than the first pass had). Two author given names on manuscript/reference.bib line 2290 are fabricated expansions of the initials that the 1971 Science byline actually carries ("M. Harada", "C. Isersky"):

1. `Isersky, Charles` -> must be `Isersky, Chaviva`.
2. `Harada, Masahiro` -> must be `Harada, Minoru`.

Corrected author line:
  author  = {Glenner, George G. and Terry, William D. and Harada, Minoru and Isersky, Chaviva and Page, David L.},

Evidence (independent of the first pass, and stronger than it):
- OpenAlex disambiguates THIS paper's own byline at the entity level: "M. Harada" -> Minoru Harada (author A5082446279) and "C. Isersky" -> Chaviva Isersky (author A5059612179). Those exact same two author IDs carry the fully spelled-out publisher metadata "Minoru Harada" / "Chaviva Isersky" on


#### `Maji2009-science` - SUPPORTS

**Claim it is cited for:** Functional amyloids mediate peptide-hormone storage ...


**Supporting passage in the cited paper:**

> I retrieved the published abstract from PubMed (PMID 19541956) and the complete body text from PMC2865899, and scanned both myself. Support for the assigned sub-clause -- peptide-hormone storage -- is the paper's title result and is unambiguous. (1) ABSTRACT, verbatim from PubMed: "We found that peptide and protein hormones in secretory granules of the endocrine system are stored in an amyloid-like cross-beta-sheet-rich conformation. Thus, functional amyloids in the pituitary and other organs can contribute to normal cell and tissue physiology." (2) RESULTS, Fig. 3D, X-ray fiber diffraction on purified membrane-depleted secretory granules -- verbatim from my own scrape of PMC2865899: "The major reflections observed were at 4.7 A interpreted as the spacing between strands in a beta-sheet and a diffuse reflection at ~10 A interpreted as the spacing between beta-sheets (Fig. 3D). These reflections are typically observed for amyloid-like fibrils." This is directly parallel to the manuscript's own WAXS framing. (3) RESULTS, Fig. 1, on the hormone-storage mechanism itself, verbatim: "...we hypothesized that ACTH might need the amyloid-forming beta-endorphin as an aggregation partner for storage in secretory granules. A 1:1 ACTH-beta-endorphin mixture in the presence of heparin formed amyloid fibrils..." Corroborating evidence in the paper: amyloid-specific antibody, Thio-T, Congo red binding and CR birefringence on purified AtT20 and rat pituitary granules (Fig. 3A-C, E-G), plus Thio-S/hormone colocalization in mouse pituitary tissue (Fig. 4). NEGATIVE RESULT I CONFIRMED MYSELF: a keyword scan of the full PMC body text returns 0 hits for "biomaterial", "mechanic", "stiff", "elastic", "robust", "hydrogel", and 0 for "extracellular matrix". Nothing in this paper speaks to the third sub-clause.


**Bibliographic corrections:** Two author given names in the .bib are factually wrong; I confirmed both independently against Crossref (10.1126/science.1173155) and the PMC2865899 byline. (1) `Perrin, Michele H.` -> `Perrin, Marilyn H.` (PMC byline: "Marilyn H Perrin, The Clayton Foundation Laboratories for Peptide Biology, The Salk Institute"). (2) `Jessberger, Sibylle` -> `Jessberger, Sebastian` (PMC byline: "Sebastian Jessberger, Institute of Cell Biology, Department of Biology, ETH Zurich"). (3) Cosmetic: `Nilsson, Karl P. R.` -> published byline is `K. Peter R. Nilsson`; "Karl" is not his published given name. IMPORTANT SEVERITY CORRECTION the first pass did not make: I compiled/inspected build/main.bbl and this manuscript uses achemso (nalefd), which prints initials only. The rendered entry reads "Maji,~S.~K.; Perrin,~M.~H.; Sawaya,~M.~R.; Jessberger,~S.; ...; Nilsson,~K. P.~R.; ...". Michele and Marilyn both in


#### `Masters1985-pnas` - SUPPORTS

**Claim it is cited for:** Pathological amyloid fibrils ... associated with Alzheimer's disease


**Supporting passage in the cited paper:**

> Two passages carry the claim. (1) INTRODUCTION, first sentence, p. 4245, left column - the association with Alzheimer's disease: "There are several closely associated morphologic changes in the brains of individuals with Alzheimer disease (AD): neurofibrillary tangles (NFTs) within neurons; plaques consisting of various proportions of amyloid cores (APCs) surrounded by neuritic degeneration; a variable degree of congophilic angiopathy (ACA); and widespread neuronal loss and gliosis in areas affected by NFTs and plaque formation." (2) RESULTS, subsection "Isolation and Purification of APCs", p. 4246, with Fig. 1c,d - the fibrillar nature: "There are two types of APC seen in the isolated state: the predominant form is the dense spherical core with a Maltese cross in polarized light (Fig. 1a), and a smaller population is composed of amorphous forms (Fig. 1b). By electron microscopy, the preparations of APC have a purity greater than 90% (Fig. 1c), and both the dense spherical cores and the amorphous cores are composed of interlacing bundles of amyloid fibrils (Fig. 1d)." The Fig. 1 caption confirms the classical amyloid identification: "Morphology of amyloid deposits: isolated and purified spherical plaque cores (a) and the amorphous variety (b) seen by polarization microscopy with Congo red ... (d) ... composed of interwoven masses of amyloid filaments." Protein identity, establishing that the AD fibril protein is a distinct one, is in the ABSTRACT ("The protein consists of multimeric aggregates of a polypeptide of about 40 residues (4 kDa)") and in the DISCUSSION, p. 4249: "The APC protein sequences reported here are not homologous to known protein sequences ... nor are they related to the sequences reported for the scrapie-associated protein." Quotes are verbatim from tesseract OCR of the PMC page scans (pnas00352-0304/0305/0307/0308), manually checked against the PMC/Crossref abstract; minor OCR artifacts in subscripts (NH2, A4) were corrected, wording was not.


#### `Palmadottir2025-bpr` - SUPPORTS

**Claim it is cited for:** The relevant materials distinction is ... whether it is effectively irreversible or can be switched by environmental conditions.


**Supporting passage in the cited paper:**

> Support is distributed across four places in the review; I read the full open-access text (Europe PMC full-text XML for PMC11836874), so all quotes below are verbatim.
> 
> (1) The "effectively irreversible" half of the claim — Sec. IV, titled "PRACTICAL IRREVERSIBILITY IN CLOSED SYSTEMS DUE TO EXPERIMENTALLY INACCESSIBLE TIME SCALES", first paragraph: "In some cases, the reverse reaction is too slow to be monitored within a practical laboratory time-frame... This has led to the interpretation in some literature as the process being irreversible although a correct interpretation would rather be that the rate constant for fibril dissociation is very low." Secs. VI and VII extend the same point ("PRACTICAL IRREVERSIBILITY IN SYSTEMS WITH CONTINUOUS MONOMER PRODUCTION"; "...WITH MONOMER CHEMICAL MODIFICATIONS").
> 
> (2) The "switched by environmental conditions" half — Sec. V, "APPARENT REVERSIBILITY UPON SYSTEM CHANGE", opening sentence: "The reversibility of amyloid formation is apparent under system change, for example, when the temperature, pH, ionic strength, pressure, denaturant concentration or solvent is altered, or upon the addition of chaperones." (citation markers removed). Fig. 6 of that section shows CD-monitored dissolution of alpha-synuclein fibrils on raising pH from 6 to 7.3, i.e. loss of the beta-sheet CD signature into random coil. Sec. VIII, "REGULATION OF FUNCTIONAL AMYLOIDS IN IN VIVO SYSTEMS BY ENVIRONMENTAL CHANGES", makes the same point biologically: functional amyloids show "growth and dissociation regulated by various factors, such as temperature, pH, nutrient availability, and developmental stage of the organism."
> 
> (3) That this dichotomy is specifically the *materials* distinction — Sec. X, "UTILIZING AMYLOID FIBRIL REVERSIBILITY FOR FUNCTIONAL BIOMATERIALS": "The fact that charged fibrils readily form hydrogels may be utilized within materials science, nanotechnology, and drug delivery, where the reversibility would serve as the key for tuning a system and/or changing the materials' properties (e.g., forming or dissolving the fibrils and the hydrogel) upon systems change." This is the single closest sentence to the manuscript's claim.
> 
> (4) The two halves are joined explicitly in Sec. XII, CONCLUSIONS: "The dissociation may, however, be so slow that the process can be viewed as practically irreversible over a given time frame or because the system at hand contains a source for constant monomer addition, leading to net buildup of amyloid even if both growth and dissociation occur at all times. Reversibility is also a key to the performance of fibril-based biomaterials and functional amyloid." The abstract carries the same two statements.
> 
> Also relevant to the fact that the distinction is a property of the assembly rather than of beta-sheet formation per se — Sec. V, closing sentences: "The different morphologies showed different heat tolerance and either dissolved, showing clear signs of reversibility upon system change or did not dissolve... Different folds will naturally possess individual stabilities and barriers and be differentially susceptible to dissolution and, thus, may or may not display reversibility over a practical laboratory time frame."


#### `Pan1993-pnas` - SUPPORTS

**Claim it is cited for:** ... prion diseases


**Supporting passage in the cited paper:**

> Three passages in the cited paper support the claim; all quoted verbatim from pdftotext extraction of the Europe PMC scan of PMC47901.
> 
> (1) DISCUSSION, p. 10965, right column -- the decisive one, and the same sentence the first pass cited AGAINST the claim: 'Although PrP amyloid plaques are diagnostic of prion diseases when present, they are often absent from both humans and animals with such diseases.' The main clause is an affirmative, and unusually strong, association statement: PrP amyloid plaques are DIAGNOSTIC of prion diseases. 'Diagnostic of' is a stronger association than the manuscript's 'associated with'.
> 
> (2) INTRODUCTION, p. 10962, right column: 'In the brains of some, but not all, animals and humans that have died of prion diseases, amyloid plaques are found which contain PrP, as determined by immunostaining and Edman protein sequencing studies (12-14).' Preceded by: 'The protease-resistant core of PrPSc designated PrP 27-30 polymerizes into rod-shaped structures which are indistinguishable from many purified amyloids both ultrastructurally and tinctorially (11).'
> 
> (3) ORIGINAL DATA -- Results, p. 10964, and Fig. 4C, restated in the Abstract: 'In contrast, polymers of PrP 27-30 were visualized as rod-shaped amyloid (Fig. 4C).' Fig. 4 legend: 'Electron micrographs of negatively stained and immunogold-labeled PrPs. (A) PrPC. (B) PrPSc. (C) Prion rods composed of PrP 27-30, negatively stained. (Bar = 100 nm.)' Plus the paper's central original result, from the Abstract: 'the beta-sheet content of PrPSc was 43% and the alpha-helix 30% as measured by FTIR', versus PrPC at 42% helix / 3% sheet.
> 
> (4) The paper explicitly places PrP in the same class as the Alzheimer amyloid precursor -- the exact grouping the manuscript's sentence makes. Discussion, p. 10966: 'In Alzheimer disease, amyloid plaques contain [beta]A4 peptide which is derived from amyloid precursor protein ([beta]-APP) that seems to undergo a similar structural transition.'


#### `Sawaya2007-nature` - SUPPORTS

**Claim it is cited for:** Many amyloid fibrils are kinetically persistent because extended hydrogen-bonded sheets are reinforced by tightly packed steric-zipper interfaces.


**Supporting passage in the cited paper:**

> I obtained the COMPLETE main text (the first pass could not — see reasoning). Four passages, all verbatim:
> 
> (1) STRUCTURAL MECHANISM — "Peptide microcrystals and fibrils", first paragraph. This is the manuscript's sentence almost word for word: "the cross-β spine consists of a pair of β-sheets; each sheet is formed from extended strands of the segment, hydrogen-bonding up and down the sheet to identical molecules, all perpendicular to the axis of the fibril. Two sheets mate tightly at a completely dry interface. At this interface, the residue side chains intermesh with close complementarity, in what we term a steric zipper." The manuscript's "extended", "hydrogen-bonded", "tightly packed", "steric-zipper" each track the paper's own wording.
> 
> (2) GENERALITY ("Many amyloid fibrils") — "Eight classes of steric zippers", first paragraph: "the structures suggest that dry, steric-zipper interfaces between β-sheets are a general principle of protein complementation in amyloid structures. Other examples of extended protein or peptide chains forming such steric-zipper interfaces between protein chains are essentially absent from the PDB ... and are rare in the CSD ..., supporting the idea that these interfaces are the defining molecular property of the amyloid state." Discussion, numbered observation (2): "The fundamental unit of amyloid-like fibrils is a steric zipper, formed by two tightly interdigitated β-sheets." Scope: ~30 segments from 14 proteins (Aβ, tau, PrP, insulin, IAPP, lysozyme, myoglobin, α-synuclein, β2-microglobulin); "in this class of proteins, we have identified one or more such segments in every protein we examined."
> 
> (3) THE STRUCTURE→PERSISTENCE CAUSAL LINK, STATED BY THE PAPER ITSELF — Figure 3 legend, final sentence. This is the passage the first pass declared did not exist: "The value of the shape complementarity parameter46, SC, for GGVVIA (SC = 0.92) is the largest value we have found for any protein interface, consistent with the higher toxicity and lower solubility of amyloid-β(1-42) than (1–40)." The authors are explicitly inferring reduced propensity to dissolve FROM tightness of zipper packing — exactly the manuscript's "because".
> 
> (4) THE PAPER REASONING FROM THE ZIPPER TO KINETICS — Discussion, numbered observation (3): "Recruitment of monomers into pre-formed fibrils is expected to be more rapid than nucleation, because recruitment requires only one molecule at a time to unmask its fibril-forming sequence, but formation of the steric-zipper nucleus requires several molecules to unmask their zipper-forming segments simultaneously. That is, the common feature of all these structures—the dry steric zipper—is compatible with slow fibril nucleation and faster fibril growth, the commonly observed kinetic characteristics of fibril formation31." (Ref. 31 is Harper & Lansbury, Annu. Rev. Biochem. 1997, "...the time-dependent solubility of amyloid proteins.")
> 
> Supporting quantitation, Supplementary Table 3 (I re-fetched and re-extracted the SI PDF independently): buried area in the dry interface 87–157 Å² per strand, 1.6–4.1 Å² per atom; dry-interface sheet-to-sheet distances 8.0–10.3 Å for the tabulated zippers; SC computed by "Lawrence and Colman's shape complementarity index" (SI footnote b, verbatim). "Tightly packed" is quantitative here, not rhetorical.
> 
> Also relevant, "Fibrils related to microcrystals": "...these segments would show enhanced protection from proton exchange... Other studies have found that three of our segments lie in regions of fibrils that show protection to proton exchange" — citing ref. 20, Kheterpal et al., "Aβ amyloid fibrils possess a core structure highly resistant to hydrogen exchange." H/D-exchange protection is a directly kinetic measure of how persistently the hydrogen-bonded core resists opening.


#### `Serpell2000-pnas` - SUPPORTS

**Claim it is cited for:** ... Parkinson's disease


**Supporting passage in the cited paper:**

> Support comes from two places in the paper, and the claim needs both halves.
> 
> (1) The Parkinson's-disease association — Abstract, first sentence (verbatim): "Filamentous inclusions made of α-synuclein constitute the defining neuropathological characteristic of Parkinson's disease, dementia with Lewy bodies, and multiple system atrophy." Reinforced in the Introduction, first paragraph (verbatim): "Parkinson's disease (PD) is the most common neurodegenerative movement disorder. Neuropathologically, it is defined by nerve cell loss in the substantia nigra and other brain regions and the presence there of Lewy bodies and Lewy neurites" and "α-synuclein is the major component of the abnormal filaments of Lewy bodies and Lewy neurites in idiopathic PD and DLB."
> 
> (2) The amyloid identity of those fibrils — Abstract, final sentence (verbatim): "X-ray diffraction and electron diffraction of the α-synuclein assemblies showed a cross-β conformation characteristic of amyloid." The underlying data are Fig. 3 (X-ray patterns from wild-type, A30P, A53T, and 1-120 α-synuclein filaments, all showing the ≈0.47 nm reflection) and Fig. 4 (electron diffraction of a partially oriented raft of wild-type filaments, 0.47 nm arc with a 0.23 nm second order on the meridian). Discussion, X-ray paragraph (verbatim): "X-ray diffraction of α-synuclein assemblies gave a cross-β pattern consisting of a 0.47-nm meridional reflection and in some samples a 1.0- to 1.1-nm equatorial reflection" and, closing that paragraph, "The fiber diffraction patterns clearly indicate that the α-synuclein filaments have a β-sheet conformation that is similar to the structure described for amyloid fibrils."
> 
> All quotes were taken verbatim from the open-access full text at https://pmc.ncbi.nlm.nih.gov/articles/PMC18329/ (Greek letters as rendered there); the abstract was independently cross-checked against the Europe PMC core record.


### Introduction paragraph 2


#### `Amiram2011-la` - SUPPORTS

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> Abstract (PMID 21258353), verbatim: "Robust high-throughput synthesis methods are needed to expand the repertoire of repetitive protein-polymers for different applications. To address this need, we developed a new method, overlap extension rolling circle amplification (OERCA), for the highly parallel synthesis of genes encoding repetitive protein-polymers... we synthesized remarkably large genes without recursive ligation. OERCA also enabled us to discover 'smart' biopolymers that exhibit fully reversible thermally responsive behaviour. This powerful strategy generates libraries of repetitive genes over a wide and tunable range of molecular weights in a 'one-pot' parallel format."


#### `Glassman2015-cw` - SUPPORTS

**Claim it is cited for:** In particular, alanine-containing ELPs can form stiff, thermoresponsive gels through arrested phase separation.


**Supporting passage in the cited paper:**

> Abstract (via Europe PMC core record, PMID 26545151). The paper reports that moderately concentrated solutions of (XPAVG)n - "where X consists of" 20% or 60% valine with the balance isoleucine - form gels rather than macrophase-separating when warmed past their inverse transition temperature. At 20 wt % in water, shear moduli reach roughly 0.1-1 MPa, with slowest stress relaxation exceeding 10^3 s. Small-angle scattering is attributed to "an arrested spinodal decomposition mechanism", modeled as a disordered bicontinuous network. The title itself is a near-verbatim statement of the manuscript claim: "Arrested Phase Separation of Elastin-like Polypeptide Solutions Yields Stiff, Thermoresponsive Gels."


#### `Glassman2016-pq` - SUPPORTS

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> Read from the free full text at PMC4752000. SEQUENCE: the gene encodes "50 repeats of the XPAVG pentapeptide, where X consists of isoleucine or valine in a 3:2 ratio"; the polymer used is "([I0.6V0.4]PAVG)50" and its chain-extended form "C-([I0.6V0.4]PAVG)50-C", with the chain tabulated as "[(IPAVGVPAVG)2(IPAVG)]10" - literal VPAVG units. MECHANISM: "The arrested phase separation of ELPs has been shown to yield remarkably stiff, biocontinuous, nanostructured networks"; above roughly 15 wt %, "the typical process of domain coalescence and macroscopic phase separation is arrested, forming a semiperiodic nanoscale network of a dense polypeptide phase"; and "polypeptide densification during phase separation arrests when a critical concentration is reached." STIFFNESS: "moduli ranging from 5 kPa to over 1 MPa over a concentration range of 5-30 wt %"; gels that "rapidly convert into stiff (G' ~ 1 MPa) and robust materials"; at 30 wt % the unextended solutions become "7 MPa gels from 0 to 37 C"; tensile Young's modulus "1.6 +/- 0.4 MPa". THERMORESPONSIVE: gelation is driven by warming through the inverse transition temperature.


#### `Hughes2018-science` - SUPPORTS

**Claim it is cited for:** ... backbone hydrogen bonding can stabilize locally ordered assemblies, whereas hydrated or weakly packed interfaces can remain dynamic or reversible.


**Supporting passage in the cited paper:**

> Verified in the PMC6192703 full text. This reference carries the sentence essentially clause by clause and is the strongest of the three.
> 
> (1) Reversible beta-sheet contacts as candidate physical junctions. Conclusion: LARKS possess three properties consistent with functioning as adhesive elements in protein gels formed from LC domains -- (i) 'High aqueous solubility contributed by their high proportion of hydrophilic residues: serine, glutamine, and asparagine'; (ii) flexibility 'ensured by their high glycine content'; (iii) multiple interaction motifs per chain, giving multivalency and 'enabling them to entangle, forming networks as found in gels.' The junction metaphor is explicit: 'If steric zippers act as molecular glue, then LARKS in LCDs act as Velcro.' The physical demonstration is Fig. 2, legend titled 'Synthetic LARKS construct forms a labile hydrogel': the 26-residue triple-LARKS peptide SYSGYSGDTSYSSYGQSNGPSTGGYG 'forms a labile hydrogel when dissolved in water at 50mg/ml and left overnight at 4 degrees C', and 'The hydrogel melts upon heating the sample to 60 degrees C for two hours', with EM confirming the fibrils melted.
> 
> (2) Backbone hydrogen bonding stabilizing locally ordered assemblies. Results, Fig. 1B-F and fig. S5: the five LARKS structures share adhesive features including 'hydrogen bonds in-register to an identical segment below it' -- the in-register backbone H-bond ladder along the fibril axis. On the inter-sheet contact: the kinks 'allow close approach of the backbones, providing favorable van der Waals or hydrogen-bond interactions between the sheets.'
> 
> (3) Hydrated OR weakly packed interfaces remaining dynamic/reversible -- BOTH limbs, contrary to the first pass's caveat. Weakly packed, explicit: 'the kinks prevent sidechains from interdigitating across the beta-sheet interface so that the kinked interfaces bury smaller surface areas than found in pathogenic amyloid fibrils, and presumably have lower binding energies'; quantified by atomic solvation parameters as '567 +/- 556 cal/mol/beta-strand' for LARKS versus '1431 +/- 685 cal/mol/beta-strand for 75 steric zipper structures', giving adhesion 'of the order of thermal energy' so that sheets 'adhere only through multivalent interactions of strands'; Fig. 1 legend supplies the metrics (shape complementarity Sc and buried solvent-accessible area Ab). Hydrated/polar, explicit in the Abstract: 'the kinked sheets interact weakly through polar atoms and aromatic sidechains' -- a polar, solvent-compatible interface, set against the amyloid case they define verbatim as 'a dry interface called a steric zipper.' Conclusion: 'the interactions must be concentration dependent and may be transient.'


**Bibliographic corrections:** ONE REAL ERROR, CONFIRMED. Author 7 in the .bib is 'Chong, Lin'; the published author is Lisa Chong. I attempted to refute this and could not -- three independent authoritative sources agree on 'Lisa': (i) the Crossref publisher deposit for 10.1126/science.aan6398 lists given name 'Lisa', family 'Chong' (ORCID 0000-0001-7667-9441); (ii) the PMC author-manuscript byline (PMC6192703, NIHMS991708) prints 'Lisa Chong'; (iii) the Europe PMC core record returns firstName 'Lisa', lastName 'Chong'. Only PubMed's compressed 'Chong L' is ambiguous, and that is almost certainly what seeded the error.
FIX: in reference.bib line ~2131, change 'Chong, Lin' to 'Chong, Lisa'.
All other fields independently re-verified correct: Hughes, Michael P.; Sawaya, Michael R.; Boyer, David R.; Goldschmidt, Lukasz; Rodriguez, Jose A.; Cascio, Duilio; Chong, Lisa; Gonen, Tamir; Eisenberg, David S. -- Science, 2018, 


#### `Kato2012-cell` - PARTIAL

**Claim it is cited for:** Reversible beta-sheet-like contacts are attractive as candidate physical junctions ...


**Supporting passage in the cited paper:**

> WHAT KATO ESTABLISHES (verified in the PMC6347373 full text).
> (1) The load-bearing main clause -- reversible beta-sheet-like contacts as candidate physical junctions of a gel. Results, Fig. 5C-D: the mCherry:FUS LC and mCherry:hnRNPA2 LC hydrogels give 'prominent, circular reflections at 4.6-4.7A and 10A', which 'are prototypic of cross-beta structure' (citing Astbury 1935; Geddes 1968; Sunde & Blake 1997), with mCherry alone and the soluble unpolymerized proteins as negative controls; the authors conclude this is 'strong evidence of the presence of amyloid-like polymers as the structural basis of hydrogel architecture.' Results, 'A Fiber Polymerization Model for Hydrogel Retention', frames the gel as 'the polymeric network of the gel.'
> (2) Reversibility/dynamics, including a genuinely THERMAL depolymerization. Summary: 'the LC sequence-based polymers described here are dynamic and accommodate heterotypic polymerization.' Fig. 6: 'mCherry:FUS LC fibers were almost fully depolymerized in the no SDS control', and -- the point the first pass missed -- 'the hnRNPA2 LC fibers were fully depolymerized upon incubation at 37 degrees C in the absence of SDS' (figure legend: 'mCherry:hnRNPA2 aggregates were fully de-polymerized under all conditions'). The b-isox precipitation is likewise explicitly thermally reversible (pellet redissolves on warming to 37 C; Fig. 1D lanes 7-9). Discussion: 'the fibers reported herein are readily de-polymerized'; they 'must be fundamentally different from the prion-like, irreversible fibers described broadly throughout the literature.'
> (3) Backbone hydrogen bonding, implicitly but adequately. Kato assigns the 4.6-4.7 A reflection to cross-beta, and cross-beta is defined by the source Kato cites (Sunde & Blake 1997) as the backbone-hydrogen-bonded strand repeat along the fibril axis. Kato's own text makes the equivalence explicit in the b-isox section: the 4.7 A spacing 'match[es] the dimensions of beta-strand - beta-strand interaction (4.7A).'
> 
> THE ONE IRREDUCIBLE GAP (basis of PARTIAL).
> Hydration / weak packing. I searched the full PMC text for water, hydration, solvation, hydrophilic, buried, interface, steric zipper, and packing. There is no analysis of interfacial hydration, buried surface area, shape complementarity, or packing density anywhere in the paper; water appears only as a methods reagent ('washed with distilled water', 'dialyzed against distilled water'). No solid-state NMR of the fibers is reported despite Tycko's authorship -- the only ssNMR mentioned is cited background on Nsp1p. The authors explicitly decline the mechanism, Discussion verbatim: 'It is clear that additional biophysical experiments will be required to understand the structural underpinnings of these newly-described amyloid-like fibers, as well as the dynamics controlling both polymerization and de-polymerization.' Their treatment of lability is phenomenological (SDS/thermal sensitivity by filtration and SDD-AGE) and comparative (versus Sup35), never interfacial.


**Why not a clean SUPPORTS:** I could not refute the first pass, and I uphold PARTIAL -- but I narrow it substantially and correct two of its supporting points, so the practical recommendation changes from 'add a reference' to 'keep as is, optionally strengthen'.

What I overturn in the first pass: (a) It claimed Kato gives no support at all for backbone hydrogen bonding, on the grounds that the words 'hydrogen bond' appear only in the b-isox crystal section. That is too literal a read. Kato's assignment of 4.6-4.7 A and 10 A to cross-beta IS the backbone-H-bond assignment as every structural biologist and the cited Sunde & Blake reference use the term, and Kato's own text equates 4.7 A with 'beta-strand - beta-strand interaction'. Clause one is supported at the level an introduction requires. (b) It characterized Kato's lability as 'detergent/dilution sensitivity at 37 C, not a thermally reversible transition'. That is wrong on the record: hnRNPA2 LC fibers were 'fully depolymerized upon incubation at 37 C in the ABSENCE of SDS', and the b-isox precipitation is explicitly reversed by warming. Kato does document thermally reversible order-disorder behavior.

What survives: the second rationale clause, 'hydrated or weakly packed interfaces can remain dynamic or reversible', is a statement about interfacial structure that Kato genuinely does not make and explicitly disclaims knowing. This is not an artifact o


#### `Luo2018-nsmb` - SUPPORTS

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> Abstract, verbatim fragments confirmed via Europe PMC and NCBI: the paper is premised on the observation that the thermostable cross-beta architecture of pathological amyloid 'does not explain' the reversibility of fibrils formed by RNA-binding proteins such as FUS. The authors identify two tandem (S/G)Y(S/G) motifs in the FUS low-complexity domain that assemble into fibrils reversing with temperature and phosphorylation state, and name them 'reversible amyloid cores, or RAC1 and RAC2', solving their atomic fibril structures by microelectron and X-ray diffraction. RAC1 is reported to form 'an ordered-coil fibril spine rather than the extended beta-strand typical of amyloids', with Ser42 -- a known FUS phosphorylation site -- essential to maintaining that ordered conformation, which supplies the mechanism for phosphorylation-dependent control of assembly. Decisively for the manuscript sentence, RAC2 forms 'a labile fibril spine with a wet interface.' The title itself is 'Atomic structures of FUS LC domain segments reveal bases for reversible amyloid fibril formation.'
> 
> Mapping to the claim: this is the reference that supplies the explicitly HYDRATED limb of 'hydrated or weakly packed interfaces can remain dynamic or reversible' -- 'wet interface' is the paper's own term, tied directly to lability, at atomic resolution. It also supplies locally ordered backbone-stabilized fibril spines (RAC1's ordered coil, RAC2's beta-strand spine) that are nonetheless thermally reversible, i.e. exactly the 'reversible beta-sheet-like contacts' of the main clause, and it is the closest structural analogue to the manuscript's thermally reversible WAXS peaks.


#### `Meyer2004-zb` - SUPPORTS

**Claim it is cited for:** ELPs ... sequence, molecular weight, and LCST behavior can be programmed recombinantly.


**Supporting passage in the cited paper:**

> Abstract, sentences 2-4 (PubMed, PMID 15132671), verbatim: "The temperature of the transition varies with ELP sequence, molecular weight, and concentration. We present a single equation of three parameters that quantitatively predicts the transition temperature as a function of ELP length and concentration for an ELP of a fixed composition. This model should be useful both for the design of new ELP sequences that exhibit a desired transition temperature and for the selection of variables to trigger the phase transition of an ELP for a given application." The first sentence names the manuscript's exact triple - sequence, molecular weight, transition behavior. The paper's central result is Tt = Ttc + (k/Length)*ln(Cc/Conc). ACS full text is genuinely unreachable (HTTP 403; Unpaywall is_oa = false with zero oa_locations; OpenAlex shows only the DOI and PubMed landing pages, no repository copy; Semantic Scholar openAccessPdf CLOSED; PMC3779073 returns "not Open Access" from the OA service). I therefore established the materials and scope from three independent, high-quality sources: (1) The Crossref reference list of Meyer 2004 itself, which contains as its ref 11 "Meyer D. E. Biomacromolecules 3 367 2002, 10.1021/bm015630n" - the authors' own recursive directional ligation paper, whose PubMed abstract (PMID 11888323) states verbatim: "We report a new strategy for the synthesis of genes encoding repetitive, protein-based polymers of specified sequence, chain length, and architecture... short gene segments are seamlessly combined in tandem using recombinant DNA techniques... We used this method to synthesize three different libraries of elastin-like polypeptides (ELPs); each library encodes a unique ELP sequence with systematically varied molecular weights... Because the thermal properties of ELPs depend on their sequence and chain length, the synthesis of these polypeptides provides an example of the importance of precise control over these parameters that is afforded by RDL." (2) McDaniel, Radford & Chilkoti, Biomacromolecules 2013, 14(8), 2866-2872 (PMC3779073, retrieved via NCBI efetch), which states verbatim: "Meyer and Chilkoti developed a quantitative model describing the effects of ELP chain length and concentration on the Tt of three different ELP libraries comprised of a mixture of Val, Ala, and Gly residues at the 4th guest residue position for a range of ELP chain lengths," and later, "Meyer and Chilkoti consolidated these observations by demonstrating that the parameters were related across three ELP libraries; they observed that k varied linearly with Ttc, whereas Cc varied as a power function with Ttc." Three libraries in 2002 = the same three sequence families characterized in 2004. (3) The field's own reading, e.g. Front. Bioeng. Biotechnol. (PMC9195583): "Extensive characterization of the mechanism underlying this phase transition permits the prediction and tuning of the LCST of the ELP, primarily by varying the hydrophobicity of the X-guest residue and the molecular weight (MW) of the ELP (Meyer and Chilkoti, 2004; McDaniel et al., 2013)" - almost verbatim the manuscript's sentence, LCST wording included, with Meyer 2004 as the citation; and a 2025 review (PMC13139353): "...results reported by McDaniel and Meyer et al. using a similar mixed ELP of (VPGXG)n, where X consisted of a mixture of A and V or A, V and G, to examine the influence of guest residue composition on ELP's Tt (McDaniel et al., 2013a; Meyer and Chilkoti, 2004)."


#### `Quiroz2015-mh` - SUPPORTS

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> Abstract (PMID 26390327), verbatim: "Yet little is known about how the phase behaviour of a protein is encoded in its amino acid sequence. Here, by synthesizing intrinsically disordered, repeat proteins to test motifs that we hypothesized would encode phase behaviour, we show that the proteins can be designed to exhibit tunable lower or upper critical solution temperature (LCST and UCST, respectively) transitions in physiological solutions. We also show that mutation of key residues at the repeat level abolishes phase behaviour or encodes an orthogonal transition... These findings set the foundation for the prediction and encoding of phase behaviour at the sequence level."


#### `Saha2020-jg` - SUPPORTS

**Claim it is cited for:** ... (both clauses, as Varanko)


**Supporting passage in the cited paper:**

> Full text read from the NIH author manuscript in PMC (PMC8297442, NIHMS1719583), retrieved as JATS XML via NCBI eutils. CLAIM (a) -- sequence/MW/LCST programmed recombinantly. Section 2, "A Brief Overview of ELPs" (with Figure 2): "ELPs are more affordably synthesized by recombinant DNA methods (Figure 2), which also allow more precise control of molecular weight and architecture." Same section, on tuning Tt: "[Urry showed] the transition temperature (T t) of ELPs could be controlled by altering the hydrophobicity of the guest residue, X, in the VPGXG motif. Building on this work, McMillan et al. determined that hydrophobic groups lower T t, while hydrophilic groups increase T t. Meyer and Chilkoti later devised an equation to account for the effect of such sequence changes as well as two other important variables -- molecular weight and concentration -- on transition temperature. This equation suggests that at a fixed pH, the chain length of the ELP is inversely related to T t." And: "McDaniel proposed a quantitative model that predicts the T t of a family of ELPs based on their composition, chain length, and concentration ... provides -- as the output -- ELPs with a specific amino acid sequence and chain length based on two inputs -- the desired T t at a specified concentration." Section 6, "ELP-Hybrid Self-Assemblies": "Their recombinant expression can be precisely controlled because of the highly regulated transcription and translation processes that create monodisperse polypeptides." CLAIM (b) -- desolvation above LCST, coacervates or hierarchical assemblies. Section 1, Introduction: "This class of biopolymers undergoes a lower critical solution temperature (LCST) phase separation in aqueous solution above a transition temperature (T t), to form an inhomogeneous coacervate -- a viscous liquid phase that is immiscible in water." Section 2: "Below its T t, an ELP adopts a random coil conformation and is well solvated; as temperature increases above its T t, the solution phase separates into an insoluble, ELP-rich coacervate phase and an ELP-poor aqueous phase." Section 3, "Self-Assembly of ELP Block Copolymers" -- the word "desolvate" is used literally: "The T t values of each block are designed to be sufficiently different to allow independent desolvation of each block. At temperatures below the critical micelle temperature (CMT), the diblock ELP is soluble. Upon raising the temperature above the CMT, the hydrophobic block selectively desolvates, which turns the diblock ELP into an amphiphile and drives its self-assembly into a micelle." Section 5, "Using Protein Order and Disorder in Hierarchical Assembly" -- this is the sentence that states the manuscript's coacervate-vs-hierarchical-assembly dichotomy directly: "ELPs alone are largely constrained to the formation of micelles or microscale liquid-like coacervates with no internal architecture. To expand on the nanostructures available to elastin-derived polymers, researchers have focused on introducing order into this otherwise disordered system." Abstract, on sequence control of architecture: "The choice of building blocks determines not only the physical properties of the nanostructures, but also their self-assembly into architectures ranging from spherical micelles to elongated nanofibers." On thermal-path dependence, Section 5 (Alpha Helices, discussing Roberts et al.): "The reversibility of phase separation could be tuned because the ordered components that drive thermal hysteresis are distinct from the disordered sequences that control the initial phase separation on heating." And Section 7 on kinetic trapping: "the authors concluded that a diblock polymer with higher hydrophobic content became kinetically trapped during the experimental time scale to form a metastable, distorted aggregate."


#### `Sing2017-zq` - SUPPORTS

**Claim it is cited for:** Above the LCST, ELPs ... form coacervates or hierarchical assemblies depending on sequence, concentration, and thermal history.


**Supporting passage in the cited paper:**

> ABSTRACT of the published paper (publisher-deposited text, retrieved via OpenAlex; RSC full text returns 403 and the article is confirmed closed access with no repository copy). Four load-bearing verbatim passages:
> 
> (1) The system is ELP thermal assembly: "Dual-associative protein di- and triblock copolymers composed of sticker-decorated midblocks and micelle-forming elastin-like polypeptide (ELP) endblocks form shear-thinning, thermoresponsively reinforceable hydrogels..."
> 
> (2) The hierarchical assembled state, and concentration + sequence-architecture dependence: "These gels, which form a disordered sphere phase due to endblock aggregation under quiescent conditions with midblock domains physically crosslinked by protein associations, exhibit both viscoelastic and thixotropic signatures with relative magnitudes dependent upon gel concentration and block architecture."
> 
> (3) Concentration dependence again, and temperature/architecture-dependent structural evolution: "For both architectures, the rate of alignment increases with increasing concentration. However, the domain formation rate when increasing temperature from 35 to 50 degrees C depends on the interplay between the thermoresponsive toughening of the endblocks and the softening of the coiled-coil domains such that the rate of rearrangement decreases in the triblock while it increases in the diblock."
> 
> (4) Explicit multiple structural length scales - i.e. hierarchy: "...recovery is characterized by a concentration-dependent restoration of the micellar network over time, with two timescales observed that correspond to two different length scales of network relaxation."
> 
> Underlying data: the quiescent in situ SAXS in Figure 2 and Scheme 1(a,b), which resolve the micelle-core / disordered-sphere-packing / crosslinked-network structure directly.


#### `Sing2018-er` - PARTIAL

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> I downloaded the green-OA full text (https://www.osti.gov/servlets/purl/1434753) and extracted the complete text with pdftotext, then grepped it exhaustively rather than relying on any summary.
> 
> SUPPORTED - alanine, and "thermoresponsive gels". Abstract, sentence 4, verbatim: "By modifying the standard glycine-containing ELP sequence (XPGVG) to instead contain alanine in the third position of the repeat sequence (XPAVG), it is possible to improve the properties of the material in both shear and extension." Alanine is the paper's controlled experimental variable (A10P4A10 vs A10P4G10 vs G10P4G10), and the materials are genuinely thermoresponsive hydrogels (15% w/v, storage moduli tracked from 20 to 50 C, end-block Tt = 14.1 C for the alanine ELP vs 24.0 C for the glycine ELP).
> 
> NOT SUPPORTED - arrested phase separation. grep over the full text: "arrest" appears exactly twice, "phase separat" exactly twice, "spinodal" zero times. Both "arrest" hits are (i) one Introduction background sentence about a DIFFERENT system, explicitly footnoted to their ref 15 = Glassman 2015 - "This substitution slows the ELP dynamics such that solutions of only ELP form tough and thermoresponsive gels due to the formation of a kinetically arrested and phase-separated nanostructure.15" - and (ii) the Glassman 2015 title in the reference list. The paper's own system is different by construction: Abstract sentence 3, "These ELPs can be used as end-blocks in triblock fusion proteins with coiled-coil associating midblock domains to result in dual-associating, network-forming materials", and Scheme 1 is captioned "...the ELP End-Blocks in (b) Results in Micellization".
> 
> STRUCTURAL POINT THE FIRST PASS MISSED, which strengthens its case: Sing 2018's nanostructure is ORDERED micellar, i.e. the structural opposite of an arrested bicontinuous network. Verbatim: "all three triblock proteins can be modeled as disordered spheres prior to flow", then under shear "the micellar... expected for cubic packing of spheres" and "micelles to a hexagonally packed structure". Glassman's arrested gels are, by contrast, "a disordered bicontinuous network" from "an arrested spinodal decomposition mechanism".
> 
> NOT SUPPORTED - stiff. "MPa" appears zero times in the entire paper. All magnitudes are kPa-scale: ultimate tensile stress "127 +/- 23 kPa", plastic plateau "around 20 kPa", steady-state shear stress "48.4 kPa" (A-type) and "0.43 kPa" (G-type). Table 1 reports yield and fracture stress in Pa. The only modulus statement is relative: "the modulus of the A-type triblock is roughly 1.5 times the moduli of the triblocks containing the glycine-rich ELP end-blocks." The paper's framing throughout is toughness, extensibility and relaxation-time scales - not stiffness.


**Why not a clean SUPPORTS:** I set out to refute the first pass and could not. Working from the full text rather than an abstract, every factual claim it made checks out, and I found one additional piece of evidence (the ordered cubic/hexagonal micellar morphology) that makes the mechanism mismatch sharper than the first pass stated. The PARTIAL verdict is UPHELD, unchanged.

The strongest available refutation is that intro citation clusters are read collectively, and that Sing 2018 does contain the sentence "...kinetically arrested and phase-separated nanostructure." I reject it for two reasons. First, that sentence is a secondary citation - Sing 2018 is restating Glassman 2015's result, with the superscript 15 attached, about "solutions of only ELP". Citing Sing 2018 for arrested phase separation is citing a paper for someone else's finding that the paper itself explicitly attributes elsewhere. Second, Sing 2018 draws the very distinction the manuscript sentence erases: its own materials are dual-associating ABA triblock fusion proteins whose gelation comes from coiled-coil midblock association plus equilibrium end-block micellization into ordered lattices. That is associative telechelic network formation, a mechanistically distinct route from arrested spinodal decomposition of a single-component ELP solution. A referee who opens this reference looking for arrested phase separation will find the paper di


#### `Varanko2020-ek` - SUPPORTS

**Claim it is cited for:** ... (same clause) AND: Above the LCST, ELPs undergo hydrophobic desolvation and can form coacervates or hierarchical assemblies ...


**Supporting passage in the cited paper:**

> Full text was read via a text-extraction proxy of the publisher page (annualreviews.org blocks direct fetch; the article is bronze OA). CLAIM (a) -- sequence/MW/LCST programmed recombinantly. Abstract: "ELPs also benefit from recombinant synthesis and genetically encoded design; these enable control over the molecular weight and precise incorporation of peptides and pharmacological agents into the sequence." Introduction, final paragraph: "This phase behavior can be controlled through selection of guest residue and chain length (9)." Section "ENGINEERING ELPs > ELP Synthesis" (with Figure 2): "ELPs are recombinant, unstructured proteins that are typically expressed in Escherichia coli"; the section then walks through concatemerization, recursive directional ligation (its ref 22 = Meyer & Chilkoti 2002, "Genetically encoded synthesis of protein-based polymers with precisely specified molecular weight and sequence"), and RDL by plasmid reconstruction, noting the latter yields "a gene with uniform length." Section "ELP Characterization": "the T t of an ELP decreases as the molecular weight of its concentration in solution increases (47). The T t is also impacted by the hydrophobicity of the guest residue in its pentapeptide motif. More hydrophobic guest residues result in a lower T t, while more hydrophilic residues increase the T t (47)." (The published sentence appears to contain a typo, "of its" for "or its"; the intended meaning is unambiguous and matches its ref 47, McDaniel, Radford & Chilkoti 2013.) CLAIM (b) -- desolvation above LCST, coacervates or hierarchical assemblies vs. sequence/concentration/thermal history. Section "ELP Characterization": "ELPs display LCST phase behavior, an entropically driven phenomenon that causes the polypeptide solution to become insoluble above its T t (9, 43) (Figure 4). This behavior is due to an unfavorable entropy of mixing, as water molecules along the polypeptide chain are highly ordered (44)." Figure 4 caption: "At temperatures below its transition temperature (T t), the ELP is hydrated and fully soluble, appearing as an optically transparent solution. As the temperature is raised above the T t of the ELP, the ELP coacervates and phase separates, appearing as a cloudy solution." Section "DRUG DELIVERY > ELP Nanoparticles" (Figure 7): "the hydrophobic block dehydrates and aggregates to form a micelle core, while the hydrophilic block remains hydrated and forms the micelle corona." Section "TISSUE ENGINEERING > Fatty Acid-Modified ELPs": "FAMEs undergo a three-stage, temperature-dependent hierarchical self-assembly that also depends on the sequence of the beta-sheet-forming peptide. At temperatures below the ELP T t, the FAME forms worm-like micelles composed of a stiff core ... and a hydrated ELP corona. Above the T t, the ELP is dehydrated and the cores aggregate, increasing interactions between nanostructures and forming mesoglobules." Section "TISSUE ENGINEERING > Partially Ordered Polymers": "The porosity and mechanical stability could be tuned by adjusting the polymer concentration and the composition or percentage of helical domains. Interestingly, the POPs remained thermally responsive but demonstrated thermal hysteresis in their phase behavior. The T t upon heating of an aqueous solution of a POP is largely a function of the ELP sequence and its chain length, but once the solid fractal network that is characteristic of the POP phase transition forms at its T t, the T t at which the POP dissolves upon cooling is largely affected by the composition and number of polyalanine domains."


### Introduction paragraph 3


#### `Bahadur2019-ln` - SUPPORTS

**Claim it is cited for:** In colloidal gels, microscopic relaxation can continue to evolve after the time-averaged structure changes only weakly ... physical aging.


**Supporting passage in the cited paper:**

> Read in full from the DOE PAGES accepted manuscript (OSTI 1574926, https://www.osti.gov/servlets/purl/1574926). SAXS and XPCS were taken in the SAME run on the SAME thermal profile -- Sec. III B opens "Simultaneous SAXS and XPCS measurements were carried out by subjecting the samples to identical thermal profiles as in the rheology studies" -- so Fig. 4 and Fig. 7 span the same t_w axis and the two halves of the claim are directly comparable.
> 
> STRUCTURE CHANGES ONLY WEAKLY -- Sec. III B (Fig. 4), 2nd paragraph: "At longer wait times there is very little change in the measured S(q) in the experimental q window." Quantified later in the same section: after the lag time the Baxter-parameter decay exponent drops "from -0.3 to -0.08 for the 82 nm particles and from -0.23 to -0.06 for the 112 nm particles," and "At longer wait times, we find that the Baxter potential doesn't change significantly and attains a plateau value." A t^-0.08 decay is a precise match to the manuscript's hedge "changes only weakly" (it is weak but nonzero, which is what the manuscript says -- not "static").
> 
> MICROSCOPIC RELAXATION CONTINUES TO EVOLVE -- Sec. III C (Fig. 7), 1st paragraph: "At even longer wait times, we observe the emergence of an intermediate Delta-g2 plateau between an initial fast decay and a slower terminal decay (not captured in these data). This plateau shifts systematically to higher values of Delta-g2 indicating that the fraction of slow relaxation modes in the system is increasing with the passage of time." Quantified from the XPCS fits (Fig. 10c,d): "the fraction of fast relaxation modes captured in 'a' shows a decrease from ~0.7 to 0.3 as a function of wait time... contributing to the exponential increase in tau_2."
> 
> THE TWO JOINED -- Sec. III B: "At t_L/t_g, the rate at which Beta changes as a function of wait time drops significantly while the elastic modulus of the bulk system starts to increase." Discussion: "At the lag time, we begin to see a divergence in the characteristic timescales of the system... We also find that the parameter 'a' keeps decreasing over time indicating that the fraction of slower relaxation modes in the system is increasing." So over t_L/t_g < t_w/t_g < 1, structure decays as t^-0.08 while tau_2 rises exponentially and the slow-mode fraction climbs from 0.3 to 0.7 -- weak structural change, strong dynamical change, both measured.
> 
> Secondary support, Discussion: "Although structural changes aren't observed at large length scales at longer wait times (t_w/t_g > 1), we find that the elastic modulus continues to evolve," attributed to rearrangements at "length scales much smaller than the particle diameter... which can't be detected using this setup."


#### `Begam2021-fh` - SUPPORTS

**Claim it is cited for:** ... protein gelation


**Supporting passage in the cited paper:**

> Mapping target: clause (ii), protein gelation. The support is in the title itself: "Kinetics of Network Formation and Heterogeneous Dynamics of an Egg White Gel Revealed by Coherent X-Ray Scattering" (verified verbatim via Crossref). The paper is an XPCS study of network formation and heterogeneous dynamics during thermal gelation of an egg-white protein gel, so it carries both halves of the manuscript's predicate — heterogeneous relaxation, and protein gelation — on its own. Independent corroboration of the mapping: Girelli et al., Nat. Commun. 2025, 16, 10814 (PMC12669764) cites this paper as its reference 32 for "gelation processes32-34" in the parallel-construction sentence quoted in the entries above.


#### `Chen2023-ho` - SUPPORTS

**Claim it is cited for:** ... (same clause)


**Supporting passage in the cited paper:**

> Read in full from the DOE PAGES accepted manuscript (OSTI 2404403, https://www.osti.gov/servlets/purl/2404403, dated 9 Dec 2022), including figure captions and the complete reference list.
> 
> AGING WITH EVOLVING MICROSCOPIC RELAXATION -- Abstract: "the suspensions form gels that undergo aging characterized by a steadily increasing elastic shear modulus and slowing, increasingly constrained microscopic dynamics." Sec. III A: "this increase persists indefinitely as the gel ages, with G' eventually growing as a weak power-law, G' proportional to t^alpha, with alpha ~ 0.4." And the XPCS observable itself: "An intermediate plateau appears at tau > 0.01 s with a plateau value that rises steadily, which indicates an increasing localization of the particles that has a close connection to the emergence of elasticity."
> 
> STRUCTURE SETTLES WHILE MICROSCOPIC RELAXATION KEEPS EVOLVING -- THE PASSAGE THE FIRST PASS MISSED. Sec. III B 2, describing Fig. 7: "The convergence time increases with increasing wave vector until about q = 0.065 nm^-1... and is almost constant at higher wave vectors... This increase in tau_X(q) with increasing q at low q indicates that the changes in the gel following the decrease in attraction strength proceed at a rate that is length-scale dependent, with LARGER-SCALE MICROSTRUCTURAL FEATURES EVOLVING MORE QUICKLY THAN SMALLER-SCALE FEATURES." The same paragraph identifies the slow, high-q observable as microscopic relaxation: "At high q, g2(q,tau) at tau = 0.1 s is in the intermediate plateau whose value reflects the localization length of the particles, which in turn can be related to the elastic modulus." That is a measured, quantitative separation of timescales in which the low-q time-averaged structure -- the same mesoscale window the manuscript's SAXS and SA-XPCS probe -- reaches its final state first while microscopic relaxation continues to evolve.
> 
> SECOND, INDEPENDENT STATEMENT OF THE STRUCTURE-QUIET HALF -- Sec. III A: "Notably, these changes in microstructure are not observed during gelation at phi = 0.43, which is close to the crossover between the gel and attractive glass states[32]." (Ref. 32 = Guo, Ramakrishnan, Harden, Leheny, J. Chem. Phys. 135, 154903 (2011) -- I confirmed this from Chen's reference list.)
> 
> TERMINOLOGY "PHYSICAL AGING" -- licensed directly. Sec. III B 1: "The observed non-monotonic behavior is instead an example of the Kovacs effect seen in aging glasses[21,22,37-41]." Reference 37 is "L. C. E. Struik, Physical Aging in Amorphous Polymers and Other Materials (Elsevier, Amsterdam, 1978)" and reference 38 is "R. L. Leheny and S. R. Nagel, 'Frequency-domain study of physical aging in a simple liquid,' Phys. Rev. B 57, 5154-5162 (1998)."
> 
> MEMORY NOT VISIBLE IN THE OBSERVABLE STRUCTURE -- Conclusion: "This reversion indicates that gels formed at even modestly different strengths of attraction have distinct and mutually incompatible microstructures despite their strongly similar properties."


#### `Chushkin2022-bw` - SUPPORTS

**Claim it is cited for:** ... cage relaxation in crowded protein solutions


**Supporting passage in the cited paper:**

> Mapping target: clause (iii), cage relaxation in crowded protein solutions. I read the complete full text (arXiv:2203.12695v1, downloaded and text-extracted; its abstract matches the published abstract word for word). Abstract, verbatim: "Here we apply an experimental design and an analysis strategy that allow us to successfully use XPCS experiments in order to measure collective long-time cage relaxation in highly crowded solutions of the eye lens protein alpha-crystallin close to and beyond dynamical arrest." The paper MEASURES stretched-exponential relaxation in the crowded near-arrest state at low dose rate. Verbatim, on phi = 0.55: "The sample at a volume fraction of phi = 0.55 is close to the glass transition but still fluid, and should thus show a stretched relaxation with beta < 1 [22, 23]. Indeed, both <tau_r> and beta show consistent and stable values below 10 kGy for the smallest two dose rates around 10 kGy/s." On phi = 0.58: "While the measurement at the lowest dose rate returns a reasonable beta < 1...". Fig. 3 caption, verbatim: "the exponent parameter beta evidences a transition from stretched exponential decay in the weakly affected fluid state to compressed decay signatures for the strongly affected state." Refs [22] and [23] invoked for beta < 1 are Bartsch, Antonietti, Schupp & Sillescu, J. Chem. Phys. 97, 3950 (1992) and Sciortino & Tartaglia, Adv. Phys. 54, 471 (2005) — the colloidal-glass literature in which beta < 1 denotes a distribution of relaxation times. Fig. 4 caption, verbatim: "The XPCS data extracted from the intrinsic dynamics in this study (red squares) are in perfect agreement with earlier measurements using rheology (data from [3]) and microrheology (data from [18])." INDEPENDENT CORROBORATION OF THE MAPPING: Girelli et al., Nat. Commun. 2025, 16, 10814 (PMC12669764) cites this exact paper as its reference 28 for "cage relaxation27,28" in the intro sentence quoted in the Girelli entry above, and again later: "Cage effects have been probed also in biological systems28,71".


#### `Girelli2021-zz` - SUPPORTS

**Claim it is cited for:** Heterogeneous relaxation resolved during protein liquid-liquid phase separation


**Supporting passage in the cited paper:**

> Mapping target: clause (i), protein liquid-liquid phase separation. Title and abstract are the primary support. Abstract (verified verbatim via Crossref/Semantic Scholar record for the DOI): "Using x-ray photon correlation spectroscopy, we determine the LLPS dynamics of a model protein solution upon low temperature quenches and find distinctly different dynamical regimes. We observe that the early stage LLPS is driven by the curvature of the free energy and speeds up upon increasing quench depth. In contrast, the late stage dynamics slows down with increasing quench depth, fingerprinting a nearby glass transition. The dynamics observed shows a ballistic type of motion, implying that viscoelasticity plays an important role during LLPS." Body content (Eq. 2 = sum of two KWW terms with separate tau1/tau2, A1/A2, gamma1/gamma2; Fig. 3(a) tau1 and tau2 vs waiting time; Fig. 3(b) nonergodicity parameter f = A2/(A1+A2) rising sharply near tw = 30 s; q-dependent exponents falling from gamma2 = 2 at low q to gamma2 = 1 at high q) is as reported by the first pass from the published PDF; I could not re-open the APS PDF myself (link.aps.org returns 403 and the DESY repository record 462418 serves no file to anonymous clients), but nothing in that reading conflicts with the abstract and I accept it as accurate. INDEPENDENT CORROBORATION OF THE MAPPING: Girelli, A. et al. Nat. Commun. 2025, 16, 10814 (open access, PMC12669764) — same senior authors — writes in its introduction, verbatim: "XPCS measurements have been demonstrated for measuring protein dynamics, both indirectly with micro-rheology25, and directly when related to Brownian diffusion26, cage relaxation27,28, liquid-liquid phase separation dynamics29-31, gelation processes32-34, and nanoscale fluctuations in hydration water35." Reference 29 there is this exact paper (Girelli et al., PRL 126, 138004 (2021)), cited for "liquid-liquid phase separation dynamics".
