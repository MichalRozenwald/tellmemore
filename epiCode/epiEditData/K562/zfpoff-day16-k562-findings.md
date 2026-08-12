# ZFPoff K562 sorted (Day 16) — CD55 methylation analysis (epicme.bio site page)

**Site page:** `docs/k562.html` (bilingual RU/EN, dark Okabe–Ito theme via `docs/assets/theme.css`). Nav entry added in `docs/index.html`. `docs/day6.html` stays reserved for the T-cell CRISPRoff-vs-unedited dataset — a *different* dataset, used here only for cross-validation.

k562.html is full-depth: per-CpG profile chart + full 137-row table, single-molecule heatmaps (50 molecules × 3 groups), per-molecule methylation-burden histograms, full ML grid (3 pairwise comparisons: LOW vs HIGH, HIGH vs BASE, LOW vs BASE), per-CpG significance (Manhattan-style) plot, top-15 lollipop chart + full top-20 ranked table, replicate concordance, TSS/KLF4 genomic mapping, T-cell cross-validation scatter, caveats, references. Both languages fully written (not machine-translated stubs).

## Data — corrected timeline (important)

The three K562 groups are **not** all from the same timepoint or procedure:

- **CD55-LOW** and **CD55-HIGH**: Day 16 post-editing, real ZFPoff, FACS-sorted on CD55 surface expression (flow). LOW: 2 reps, 4,333 molecules. HIGH: 2 reps, 2,868 molecules.
- **Baseline**: K562 at **Day 2** of passaging, post-**mock** ("fake"/control) editing — not real ZFPoff, not the same day as LOW/HIGH. 1 rep, 206 molecules. Two uploaded versions of this file differ ~2x in methylated-call count (14.8% vs 7.4% mean methylation, same 206×137 shape); the page uses the 14.8% version (matched 3 of 4 uploads) — unconfirmed which is correct.

**Consequence:** any comparison involving baseline (LOW-vs-baseline, HIGH-vs-baseline) confounds "edited vs not" with a large timepoint/passaging/procedure difference (Day 2 mock vs Day 16 real ZFPoff). **LOW vs HIGH is the only clean same-day, same-procedure comparison** (both Day 16, both real ZFPoff, differ only by sort gate) — it's the one the page (and this doc) leans on for the headline finding.

- Single-molecule CpG methylation matrices: reads × 137 padded CpG positions, 1 = methylated call, NaN = not called methylated.
- Genomic coordinates (GRCh38, chr1:207,318,092–207,324,450) obtained from `20260716_mCG_COORDINATED_FRACTIONS_D6_CROff_vs_Unedited_Tcells_0995thresh_filtered.csv` — **this coordinates file is from the T-cell dataset, not K562** (confirmed by user). It supplied an independent T-cell CRISPRoff-vs-unedited comparison for the same 137 CpGs, used here only for cross-cell-type validation (the site's existing `docs/day6.html` T-cell page).

## Key findings

- **CD55-HIGH cells are hypermethylated relative to CD55-LOW** at this locus: 39.9% mean CpG methylation (HIGH, Day 16) vs 12.9% (LOW, Day 16). Baseline (Day 2, mock) sits at 14.8%, close to LOW — but given the timepoint/procedure mismatch, this proximity is not strong evidence either way. Counter-intuitive vs. the canonical promoter-methylation-silences-expression model, which would predict LOW to be hypermethylated.
- ML classifiers (logistic regression + random forest, 5-fold CV) separate CD55-LOW vs CD55-HIGH molecules with AUC ≈ 0.92–0.93 — this is the load-bearing, timepoint-matched comparison. HIGH-vs-baseline (AUC ≈ 0.95–0.96) and LOW-vs-baseline (AUC ≈ 0.68–0.71) are reported but should be read as edited-vs-mock-at-a-different-timepoint, not a clean edited-vs-unedited contrast.
- Per-position Fisher's exact test + BH-FDR: 119/137 CpGs differ LOW vs HIGH (q<0.05); 125/137 HIGH vs baseline; only 18/137 LOW vs baseline.
- The differential signal (LOW vs HIGH) clusters at CpGs #49–122 (−374 to +757 bp relative to the CD55 TSS, chr1:207,321,678), peaking at CpGs #78–83 (+55 to +84 bp, just downstream of TSS). CpGs #35, #54, #65 sit within ~15 bp of the three published KLF4 CACCC-box sites (Site C −667/−652, Site B −310/−302, Site A −106/−95; Zhang et al., PMC10884306) — confirms this amplicon covers the CD55 core promoter/TSS region.
- **Cross-cell-type validation:** the same per-CpG delta pattern (HIGH−LOW in K562, Day 16) correlates r=0.777 with the independent T-cell dataset's delta (CRISPRoff − unedited, Day 6) — the hypermethylation-with-HIGH-expression direction reproduces outside K562/sorting and outside this specific timepoint.

## Open todos

- Confirm which baseline (RVS118) upload is correct (14.8% vs 7.4% global methylation, same 206×137 shape).
- Investigate why CD55-HIGH (not LOW) carries the methylation — framed as an open question specifically about the LOW-vs-HIGH (Day 16, same-procedure) result, independent of the baseline dating issue.
- Consider a true same-day K562 unedited control (if one exists) to properly isolate the "was this cell edited at all" effect, since current baseline is confounded by both time and mock-vs-real procedure.
- `about.html`, `about-cd55.html`, `day6.html`, `day35.html`, `day6VSday35.html` haven't been shared yet — needed for any further cross-page alignment work.

## Sources
- KLF4/CD55 promoter mechanism: https://pmc.ncbi.nlm.nih.gov/articles/PMC10884306/
- CD55 genetic/epigenetic regulation, colon cancer (TFCP2, NF-κB, miR-27a-3p): https://pmc.ncbi.nlm.nih.gov/articles/PMC9889927/
- GeneCards CD55: https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD55
- NCBI Gene CD55 (GRCh38 coordinates): https://www.ncbi.nlm.nih.gov/gene/1604
