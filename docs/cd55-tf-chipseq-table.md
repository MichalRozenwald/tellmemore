# CD55 (DAF) — transcription factors and available ChIP-seq data

Reference table for epiCme. Ordered by the published CD55-regulation literature first,
then by mechanistic importance for the K562 CD55 promoter window
(chr1:207,318,092–207,324,450, GRCh38; TSS at chr1:207,321,678).

**Legend**
- **Accession in bold** = verified ENCODE experiment, checked individually.
- "All ENCODE" = filter query returning every released TF ChIP-seq for that target, all cell lines.
- "In K562 ✓" = peak confirmed in the CD55 window per UCSC ENCODE cluster tracks.

ENCODE filter URL pattern:
`https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=<TF>`

---

## Part 1 — Published CD55 regulators (order as in the source list)

### 1. NF-κB — RELA

| Field | Value |
|---|---|
| Verified ENCODE | **ENCSR000EAG** — GM12878 + TNFα 6 h — https://www.encodeproject.org/experiments/ENCSR000EAG/ |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=RELA |
| Other cell lines | Abundant across HepG2, K562, HeLa-S3, lymphoid lines |
| K562 / T-cell status | Lymphoid data exists. ENCODE K562 clusters are **unstimulated** — NF-κB is stimulus-dependent, so a missing peak is a condition mismatch, not evidence of no binding |
| Role in CD55 regulation | Promoter mutagenesis + luciferase show dependence on an NF-κB site; ChIP supports binding at the CD55 promoter in colon cancer cells |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC9889927/ |

### 2. NF-κB — NFKB1 (p50)

| Field | Value |
|---|---|
| Verified ENCODE | not individually verified |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=NFKB1 |
| Other cell lines | GM12878, K562 experiments exist |
| K562 / T-cell status | Same stimulus caveat as RELA |
| Role in CD55 regulation | p50 subunit of the same NF-κB complex implicated at the promoter |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC9889927/ |

Related: **RELB** (non-canonical NF-κB), GM12878 — **ENCSR387QUV** —
https://www.encodeproject.org/experiments/ENCSR387QUV/

### 3. TFCP2 (CP2 / LBP-1 family)

| Field | Value |
|---|---|
| Verified ENCODE | none found |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=TFCP2 |
| Other cell lines | **No ENCODE experiment located in any cell line.** Only realistic sources: GTRD, ChIP-Atlas, or the source paper's own ChIP |
| K562 / T-cell status | Not confirmed anywhere; the single genuine dead end in this list |
| Role in CD55 regulation | Knockdown, site mutation, and ChIP support direct TFCP2 binding to the CD55 promoter and regulation of CD55 expression |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC9889927/ |

### 4. KLF4 (KLF family, CACCC motifs)

| Field | Value |
|---|---|
| Verified ENCODE | none in a native context |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=KLF4 |
| Other cell lines | ENCODE has eGFP-tagged KLF4 in HEK293 only. GEO: vascular endothelial KLF4 ChIP-seq (laminar-shear-stress studies); iPSC reprogramming time-course series; mouse ESC GSE11431 |
| K562 / T-cell status | **Not profiled in K562 or T cells.** Absence carries no evidence about CD55 — KLF4 was simply never an ENCODE K562 priority target |
| Role in CD55 regulation | Three CACCC-box sites near the TSS; the leading proposed activator. Hypothesis remains motif-based and functional, not ChIP-confirmed in K562 |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC10884306/ |

### 5. CREB (via CRE sites) — CREB1

| Field | Value |
|---|---|
| Verified ENCODE | **ENCSR000BSO** — K562 — https://www.encodeproject.org/experiments/ENCSR000BSO/ |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=CREB1 |
| Other cell lines | Abundant: GM12878, HepG2, LNCaP, many more |
| K562 / T-cell status | K562 experiment exists. Whether an IDR peak falls inside the CD55 window still needs a direct bigBed check — a prior summarized-API pass did not confirm it |
| Role in CD55 regulation | CD55 promoter constructs with CREB-site mutations show reduced induction in an HCV-core context; CRE/CREB signalling also implicated alongside KLF4. Note **ATF1** (same family, same CRE) is confirmed in the window — see Part 2 |
| References | https://pmc.ncbi.nlm.nih.gov/articles/PMC3558599/ · https://pmc.ncbi.nlm.nih.gov/articles/PMC10884306/ |

### 6. SP1

| Field | Value |
|---|---|
| Verified ENCODE | **ENCSR000BKO** — K562, Myers/HAIB — https://www.encodeproject.org/experiments/ENCSR000BKO/ |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=SP1 |
| Other cell lines | Very abundant across GM12878, HepG2, HeLa-S3, HEK293 |
| K562 / T-cell status | K562 experiment confirmed. Window-level peak not yet confirmed — surprising for a CpG-island promoter, so worth re-checking from raw bigBed |
| Role in CD55 regulation | SP1 binding is functionally required for CD55 promoter activity regulation |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC3558599/ |

### 7. HIF-1 (HIF1A)

| Field | Value |
|---|---|
| Verified ENCODE | none |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=HIF1A |
| Other cell lines | **GEO GSE145157** — HIF1A ChIP-seq in HCT116, RKO, A549, H460 after 24 h hypoxia · MCF-7 HIF ChIP-seq (Blood 2011) · DLD-1 and TIG-3, 531 and 616 sites · GSE16347 HepG2 (ChIP-chip) |
| K562 / T-cell status | Not in K562 or T cells. Requires hypoxic induction, so normoxic profiling would be uninformative anyway |
| Role in CD55 regulation | Hypoxia work reports HIF-dependent induction and describes putative HIF-1 binding sites in the CD55 promoter |
| References | https://ashpublications.org/blood/article/117/23/e207/22454/High-resolution-genome-wide-mapping-of-HIF-binding · https://www.nature.com/articles/s41467-021-21687-2 · https://pmc.ncbi.nlm.nih.gov/articles/PMC3051044/ |

### 8. PU.1 (SPI1)

| Field | Value |
|---|---|
| Verified ENCODE | **ENCSR000BGQ** — GM12878 — https://www.encodeproject.org/experiments/ENCSR000BGQ/ |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=SPI1 |
| Other cell lines | Abundant in myeloid and B-lineage lines; K562 likely available |
| K562 / T-cell status | Lymphoid confirmed. PU.1 is largely silent in mature T cells, so a T-cell negative would be biologically expected |
| Role in CD55 regulation | Allele-dependent CD55 promoter effect linked to PU.1 binding; SPI1 knockdown reduces promoter activity in monocyte models |
| Reference | **TO DO — no PMID supplied in the source list; needs pinning before publication** |

---

## Part 2 — Other TFs, ordered by importance for CD55 regulation

Ranking rationale: factors confirmed in the K562 CD55 window and mechanistically
capable of explaining the LOW/HIGH methylation difference come first; general
transcriptional-machinery components that mark the region as active come last.

### 9. CTCF — top mechanistic candidate

| Field | Value |
|---|---|
| Verified ENCODE | **ENCSR000EGM** — K562, Snyder lab — https://www.encodeproject.org/experiments/ENCSR000EGM/ |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=CTCF |
| Other cell lines | Profiled in essentially every ENCODE line |
| K562 status | In K562 ✓ — peak at CpG #43–86, chr1:207,320,854–207,321,372, inside the significant cluster |
| Role in CD55 regulation | Not a literature-reported CD55 regulator, but binding is directly blocked by CpG methylation in its recognition site. Rising methylation in CD55-HIGH could displace it — a direct mechanism rather than a correlate |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC10325916/ |

### 10. NRF1 — top mechanistic candidate

| Field | Value |
|---|---|
| Verified ENCODE | not individually verified; confirmed in K562 via UCSC cluster tracks |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=NRF1 |
| Other cell lines | K562, GM12878, HepG2 and others |
| K562 status | In K562 ✓ — peak at CpG #13–43, chr1:207,320,698–207,321,177 |
| Role in CD55 regulation | The textbook methylation-sensitive factor: binding is directly competed by CpG methylation. Same displacement hypothesis as CTCF |
| Reference | https://www.nature.com/articles/nature16462 |

### 11. ATF1 — the actual ChIP support for the "CREB" claim

| Field | Value |
|---|---|
| Verified ENCODE | not individually verified; confirmed in K562 via UCSC cluster tracks |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=ATF1 |
| Other cell lines | K562, GM12878, HepG2 |
| K562 status | In K562 ✓ — within the CD55 promoter window |
| Role in CD55 regulation | CREB/ATF family member binding the same TGACGTCA CRE and heterodimerising with CREB1. The CRE module implicated by the CD55 literature **is** occupied in K562 — under a paralog's name |
| Reference | https://pmc.ncbi.nlm.nih.gov/articles/PMC3558599/ |

### 12. DNMT3B — the candidate writer

| Field | Value |
|---|---|
| Verified ENCODE | not individually verified |
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=DNMT3B |
| Other cell lines | **HepG2 only** for the peak in this window |
| K562 status | No K562 experiment. Cannot be confirmed or excluded in this line |
| Role in CD55 regulation | De novo methyltransferase — would be the enzyme establishing the LOW→HIGH methylation difference if confirmed. Currently indirect |
| Reference | https://www.nature.com/articles/nrg.2016.83 |

### 13. EP300 / p300 — activity context

| Field | Value |
|---|---|
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=EP300 |
| K562 status | In K562 ✓ |
| Role in CD55 regulation | Coactivator, not sequence-specific. Confirms the window is an active regulatory module rather than a silent region |
| Reference | — |

### 14. POLR2A (total and Ser5-P) — activity context

| Field | Value |
|---|---|
| All ENCODE | https://www.encodeproject.org/search/?type=Experiment&assay_title=TF+ChIP-seq&status=released&target.label=POLR2A |
| K562 status | In K562 ✓, both total and initiating Ser5-P forms |
| Role in CD55 regulation | Engaged and initiating Pol II at the CD55 TSS. Establishes that the methylation change sits in a transcriptionally active promoter |
| Reference | — |

---

## Open items before publication

1. **TFCP2** — check GTRD and ChIP-Atlas. If both return nothing, state the absence
   explicitly as a negative result.
2. **KLF4** — same; the endothelial and reprogramming GEO datasets may still be
   intersected with the CD55 window even though the cell type is wrong.
3. **SPI1 reference** — no PMID was supplied for the PU.1/DAF allele claim. Must be
   pinned rather than inferred.
4. **SP1 and CREB1 window-level peaks** — re-check from raw ENCODE bigBed
   (`hgdownload.soe.ucsc.edu` + `bigBedToBed`) rather than the summarized API.
5. **Cross-assembly** — ENCODE is hg38 only. liftOver hg38 → hs1 (T2T) for any peak
   set used alongside T2T coordinates.

## Aggregator databases for the remaining gaps

- GTRD — http://gtrd.biouml.org — 15,982 human ChIP-seq experiments, 1,391 TFs,
  uniformly processed, per-TF metaclusters, search-near-gene function
- ChIP-Atlas — https://chip-atlas.org/peak_browser
- ReMap 2022 — https://remap.univ-amu.fr/
- Cistrome DB — http://cistrome.org/db/
- ChIPSummitDB — http://summit.med.unideb.hu/summitdb/ — motif-to-summit distances,
  useful for distinguishing direct binding from co-IP
- UniBind — https://unibind.uio.no/ — ChIP-supported and motif-matched only
