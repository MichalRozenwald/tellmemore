# ZFPoff K562 sorted (Day 16) — CD55 methylation analysis (epicme.bio site page)

**Site page:** `docs/k562.html` (bilingual RU/EN, dark Okabe–Ito theme via `docs/assets/theme.css`). Nav entry added in `docs/index.html`. `docs/day6.html` stays reserved for the T-cell CRISPRoff-vs-unedited dataset — a *different* dataset, used here only for cross-validation.

## Data pipeline correction (important, changes the headline numbers)

The two uploaded files per group (LOW, HIGH, BASE) are **not biological/technical replicates**. They are the *same* Oxford Nanopore sequencing run, methylation-called at two different confidence thresholds for deciding whether a given C counts as methylated: **0.7 (permissive)** and **0.995 (strict)**. An earlier version of this analysis pooled (concatenated) both files per group as if they were independent replicates — this double-counts the same physical molecules and inflates N, which is statistically invalid (violates the independence assumption behind Fisher's exact test and the ML cross-validation).

**Fix:** every headline number now uses only the **0.995 (strict) threshold** file per group — the same threshold already used in the T-cell validation dataset (its coordinates file has `mc0995` columns), so the cross-dataset comparison is apples-to-apples. The 0.7 threshold file is kept only as a "threshold sensitivity" check (Pearson correlation of the per-CpG profile between 0.7 and 0.995), not pooled into the primary analysis. This also resolves the earlier "two baseline uploads, unclear which is correct" question: both are correct — they're the same 206 molecules at two different calling thresholds (14.8% at 0.7, 7.4% at 0.995).

### Corrected numbers (0.995 threshold, single file per group, no pooling)

- **N molecules:** LOW=2,163, HIGH=1,433, BASE=206
- **Global % methylation:** LOW=9.2%, HIGH=29.4%, BASE=7.4%
- **LOW vs HIGH:** 117/137 CpGs significant q<0.05; LR AUC=0.922, RF AUC=0.912
- **HIGH vs BASE:** 114/137 significant; LR AUC=0.949, RF AUC=0.938
- **LOW vs BASE:** 18/137 significant; LR AUC=0.526, RF AUC=0.524 — near-chance once the pooling artifact is removed. CD55-LOW (Day 16, low methylation) is barely distinguishable from baseline (Day 2, mock); HIGH is the outlier that gained substantial methylation.
- **T-cell cross-validation r:** 0.784
- **Threshold-sensitivity correlation** (0.995 vs 0.7 profile, same group): LOW r=0.976, HIGH r=0.906, BASE r=0.974 — the qualitative pattern is fairly robust to threshold choice.
- **Significant CpG cluster (LOW vs HIGH):** peak at CpG #78–83 (+58 to +84 bp from TSS); main contiguous significant block at #20–132 (−1,077 to +1,935 bp from TSS).

## Group labels — LOW/HIGH are methylation level, not expression level

**CD55-LOW and CD55-HIGH denote CD55 methylation level, not CD55 expression level.** Cells were flow-sorted using CD55 surface expression as the readout, but the resulting fractions are named for their methylation status: CD55-LOW = the low-methylation fraction, CD55-HIGH = the high-methylation fraction. Under the canonical model (promoter methylation → gene silencing), the low-methylation (LOW) fraction is expected to correspond to *higher* surface CD55 expression, and the high-methylation (HIGH) fraction to *lower* surface expression — i.e., HIGH being more methylated than LOW is the expected canonical direction, not a contradiction of it. This dataset does not have a direct, paired single-cell measurement of expression alongside methylation — the inverse relationship follows from the model and the sorting method, not from a direct joint readout.

## Data — timeline

The three K562 groups are **not** all from the same timepoint or procedure:

- **CD55-LOW** and **CD55-HIGH**: Day 16 post-editing, real ZFPoff, FACS-sorted on CD55 surface expression as a readout of methylation level.
- **Baseline**: K562 at **Day 2** of passaging, post-**mock** ("fake"/control) editing — not real ZFPoff, not the same day as LOW/HIGH.

**Consequence:** any comparison involving baseline confounds "edited vs not" with a large timepoint/passaging/procedure difference. **LOW vs HIGH is the only clean same-day, same-procedure comparison** — the one the page leans on for the headline finding.

- Genomic coordinates (GRCh38, chr1:207,318,092–207,324,450) obtained from `20260716_mCG_COORDINATED_FRACTIONS_D6_CROff_vs_Unedited_Tcells_0995thresh_filtered.csv` — this coordinates file is from the T-cell dataset, not K562. It supplied an independent T-cell CRISPRoff-vs-unedited comparison for the same 137 CpGs, used here only for cross-cell-type validation (the site's existing `docs/day6.html` T-cell page — "CRISPRoff" is confirmed the correct term for that dataset).

## Key findings

- CD55-HIGH (high methylation, 29.4%) vs CD55-LOW (low methylation, 9.2%) — this direction is the expected, canonical-model-consistent pattern once LOW/HIGH are correctly understood as methylation-level labels.
- LOW vs HIGH is the load-bearing, timepoint-matched comparison (AUC≈0.92). LOW vs BASE is near-chance (AUC≈0.53) once the pooling artifact is removed — LOW looks much like baseline; HIGH is the outlier.
- Differential signal (LOW vs HIGH) peaks at CpGs #78–83 (+58 to +84 bp from TSS, chr1:207,321,678), with a broader significant block at #20–132. CpGs #35, #54, #65 sit within ~15 bp of the three published KLF4 CACCC-box sites (Site C −667/−652, Site B −310/−302, Site A −106/−95; Zhang et al., PMC10884306) — confirms this amplicon covers the CD55 core promoter/TSS region.
- Cross-cell-type validation: the K562 HIGH−LOW delta correlates r=0.784 with the independent T-cell CRISPRoff−unedited delta (Day 6) — the direction reproduces outside K562/sorting and outside this specific timepoint.

## Open todos

- Get a direct, paired single-cell measurement of CD55 surface expression alongside methylation for the LOW/HIGH fractions, to directly confirm (rather than infer from the sort/model) the expected inverse relationship.
- Consider a true same-day K562 unedited control (if one exists) to properly isolate the "was this cell edited at all" effect, since baseline is confounded by both time and mock-vs-real procedure.
- `about.html`, `about-cd55.html`, `day6.html`, `day35.html`, `day6VSday35.html` haven't been shared yet — needed for any further cross-page alignment work.

## Sources
- KLF4/CD55 promoter mechanism: https://pmc.ncbi.nlm.nih.gov/articles/PMC10884306/
- CD55 genetic/epigenetic regulation, colon cancer (TFCP2, NF-κB, miR-27a-3p): https://pmc.ncbi.nlm.nih.gov/articles/PMC9889927/
- GeneCards CD55: https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD55
- NCBI Gene CD55 (GRCh38 coordinates): https://www.ncbi.nlm.nih.gov/gene/1604
