# topGO over-representation analysis — corrected pipeline

`topGO_fixed.r` replaces `topGO.r`. The original script had two defects that
between them removed the strongest biological signal from the results and made
the multiple-testing correction meaningless. Both are described below with the
evidence, followed by a paper-ready methods paragraph and the figure legends.

`topGO.r` and its outputs (`overall_clay_vs_sandy/GO_up.*`, `GO_down.*`) are kept
unmodified for comparison. New results are written to
`overall_clay_vs_sandy/topgo_corrected/`.

---

## 1. What was wrong, and what changed

### Bug 1 — FDR correction applied to the already-significant subset

```r
table1 <- filter(table, Classic < 0.05)                # keep only p < 0.05
p.adj  <- round(p.adjust(table1$Classic, "BH"), 4)     # ...then "correct" within that subset
```

Benjamini–Hochberg must be computed over **all** tested hypotheses. Applying it
to a pre-filtered set makes the correction vacuous: the largest raw p in the
subset is already < 0.05, and BH's largest adjusted value is `p_max × n/n =
p_max`. Every surviving term therefore passes "FDR < 0.05" by construction.

In the old `GO_up.csv` all 170 rows passed, and the last row's adjusted p
(0.0472) is identical to its raw p — the signature of the bug. `results.table.bh`
and `results.table.p` were always the same object.

This also explains the near-constant adjusted p-values that prompted the review:
BH assigns adjusted values as a running minimum from the bottom of the ranked
list upward, so a small `n` collapses long runs to a single number. The old file
had 31 consecutive terms at 0.0035 and 39 at 0.0339.

**Fix:** take numeric scores for every tested term from `score(result)` and run
`p.adjust()` over all of them.

### Bug 2 — p-values compared as text

`GenTable()` returns its score column as **character**, formatted for printing.
So `filter(table, Classic < 0.05)` performed a lexicographic string comparison:

```r
"0.00011" < "0.05"   # TRUE
"4.1e-09" < "0.05"   # FALSE   <- scientific notation sorts after "0"
```

Every p-value small enough to be printed in scientific notation (below ~1e-4) was
silently discarded. **The 31 most significant up-regulated terms never reached
the output file**, including `carboxylic acid metabolic process` (p = 4.1e-09),
`oxoacid metabolic process`, `small molecule metabolic process`, and the entire
glycolysis / gluconeogenesis block. The best p-value in the old `GO_up.csv` was
0.00011; the true best is 4.1e-09.

**Fix:** never parse `GenTable()`'s formatted strings — read numeric scores from
the result object and use `GenTable()` only for term names and counts.

### Other changes

| | `topGO.r` | `topGO_fixed.r` |
|---|---|---|
| `nodeSize` | unset (all 5,886 BP terms tested, including 1–2 gene terms) | `5` (2,633 terms) |
| Term names | truncated at 40 characters (`purine ribonucleoside diphosphate metabo...`) | `numChar = 1000`, full names |
| Algorithms | `classic` only | `classic` + BH, and `weight01` alongside |
| up / down | line edited by hand; output names hardcoded to `GO_down` | single loop, names derived from the gene set |
| CSV quoting | `quote = FALSE` — GO terms containing commas shifted every later column | default quoting |
| Paths | `/home/diegoj/...` (no longer exists) | `/dados04/jorge/...` |

### Effect on the results

| | old `GO_up.csv` | corrected |
|---|---|---|
| BP terms tested | 5,886 | 2,633 |
| best raw p reported | 0.00011 | **4.1e-09** |
| terms at FDR < 0.05 | 170 (100% of nominal — meaningless) | **67** |
| terms at weight01 p < 0.05 | not computed | 28 |

For the **down-regulated** set the corrected analysis returns **no significant
term at all** (0 of 2,633 at FDR < 0.05). The old `GO_down.csv`, which listed 62
"significant" terms, was entirely an artefact: that gene set contains only 5
genes, 4 of them annotated, and every term it produces rests on a single gene.

---

## 2. Methods

> Gene Ontology over-representation analysis was performed with topGO v2.x
> (Alexa & Rahnenführer) on the differentially expressed gene sets
> (|log2FC| > 1, adjusted p < 0.05; DESeq2). The gene universe comprised the
> 3,024 genes carrying at least one Biological Process annotation in the
> reference annotation, and terms with fewer than five annotated genes were
> excluded (`nodeSize = 5`), leaving 2,633 testable terms. Enrichment was
> assessed by Fisher's exact test under the `classic` algorithm, with
> Benjamini–Hochberg correction applied across all tested terms; terms with
> FDR < 0.05 were considered significant. Because `classic` scores each term
> independently of the GO graph and therefore returns redundant ancestor–
> descendant terms, results were additionally computed with the topology-aware
> `weight01` algorithm, which scores each term conditional on its descendants.
> `weight01` p-values are reported uncorrected, as the algorithm already
> accounts for the dependency between terms and further correction is not
> recommended. Figures use the perceptually uniform, colour-vision-deficiency
> safe scientific colour map *batlow* (Crameri 2018).

**References**

- Alexa A, Rahnenführer J, Lengauer T (2006) Improved scoring of functional
  groups from gene expression data by decorrelating GO graph structure.
  *Bioinformatics* 22:1600–1607.
- Alexa A, Rahnenführer J. topGO: Enrichment Analysis for Gene Ontology.
  R package.
- Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate.
  *J R Stat Soc B* 57:289–300.
- Crameri F (2018) Scientific colour maps. Zenodo, doi:10.5281/zenodo.1243862.
- Crameri F, Shephard GE, Heron PJ (2020) The misuse of colour in science
  communication. *Nat Commun* 11:5444.

---

## 3. Figure legends

### `GO_up_classic_BH.{pdf,svg,png}`

> **GO biological processes over-represented among genes up-regulated in clay
> relative to sandy soil.** Each point is a Gene Ontology Biological Process
> term; the 20 most significant of the 67 terms passing FDR < 0.05 are shown,
> ordered top to bottom by increasing adjusted p-value. Horizontal position gives
> the Benjamini–Hochberg adjusted p-value on a −log10 scale, so points further
> right are more significant. Point area encodes the number of up-regulated genes
> annotated to the term (4–24 genes) and fill colour the fold enrichment over the
> count expected by chance (3.0–44.8×). The dashed vertical line marks
> FDR = 0.05; all points lie to its right. Over-representation was tested by
> Fisher's exact test under the topGO `classic` algorithm against a universe of
> 3,024 annotated genes, and p-values were corrected across all 2,633 Biological
> Process terms tested. Because `classic` scores every term independently of the
> GO graph, ancestor and descendant terms that share the same genes appear as
> separate entries: the three leading acid-metabolism terms rest on the same 19
> genes, and the six glycolysis-related terms in the lower half of the panel on
> the same 4, all six returning an identical p-value.

### `GO_up_weight01.{pdf,svg,png}`

> **Topology-corrected GO biological processes over-represented among genes
> up-regulated in clay relative to sandy soil.** Each point is a Gene Ontology
> Biological Process term; the 20 most significant of the 28 terms with
> weight01 p < 0.05 are shown, ordered top to bottom by increasing p-value.
> Horizontal position gives the weight01 p-value on a −log10 scale, point area
> the number of up-regulated genes annotated to the term (2–8 genes), and fill
> colour the fold enrichment over the count expected by chance (7.3–44.8×). The
> dashed vertical line marks p = 0.05. The topGO `weight01` algorithm traverses
> the GO graph from the most specific terms upward and scores each term
> conditional on its descendants, so broad terms whose signal is already
> explained by a more specific child are down-weighted: 102 of the 137 terms
> nominally significant under the `classic` algorithm fall to p = 1 here. This is
> why general terms such as *carboxylic acid metabolic process* are absent while
> the specific term *canonical glycolysis* is retained, and why this panel gives
> a non-redundant view of the same response. p-values are reported without
> multiple-testing correction, as the algorithm conditions each test on the graph
> structure and thereby already accounts for the dependency between terms. Gene
> universe: 3,024 annotated genes; 2,633 Biological Process terms tested.

### `GO_down_weight01.{pdf,svg,png}`

> **Exploratory GO biological process results for genes down-regulated in clay
> relative to sandy soil; no term is statistically significant.** No Gene
> Ontology term reached significance after FDR correction in this gene set
> (0 of 2,633 Biological Process terms at FDR < 0.05). The down-regulated set
> contains only 5 genes, 4 of which carry a GO annotation. All 19 terms with an
> uncorrected weight01 p-value below 0.05 are shown, ordered top to bottom by
> increasing p-value. Horizontal position gives the weight01 p-value on a −log10
> scale, point area the number of down-regulated genes annotated to the term, and
> fill colour the fold enrichment over the count expected by chance. **Every term
> shown is supported by a single gene** (all points share the smallest size in
> the legend), so the apparent fold enrichments of 21.6–151.2× are an artefact of
> very small counts and must not be interpreted as biological enrichment. The
> dashed vertical line marks p = 0.05. This panel is presented for completeness
> only. Gene universe: 3,024 annotated genes.

### `GO_down_classic_rawp_NOT_significant.{pdf,svg,png}`

> **Exploratory GO biological process results for genes down-regulated in clay
> relative to sandy soil, uncorrected p-values; no term is statistically
> significant.** No Gene Ontology term survived FDR correction in this gene set
> (0 of 2,633 Biological Process terms at FDR < 0.05), so the 20 terms with the
> smallest **uncorrected** p-values under the topGO `classic` algorithm are shown
> instead, ordered top to bottom by increasing p-value. Horizontal position gives
> the uncorrected Fisher's exact p-value on a −log10 scale, point area the number
> of down-regulated genes annotated to the term, and fill colour the fold
> enrichment over the count expected by chance. The down-regulated set contains
> only 5 genes, 4 of which carry a GO annotation, and **every term shown is
> supported by a single gene** (all points share the smallest size in the
> legend); the apparent fold enrichments of 47–151× are an artefact of very small
> counts. The dashed vertical line marks p = 0.05. No claim of enrichment is made
> from this panel; it is shown for completeness only. Gene universe: 3,024
> annotated genes.

---

## 4. Running it

```bash
cd /dados04/jorge/bianca/rnaseq/full_run1_no_collapse/star_salmon/deseq2_qc
/home/genomics/miniconda3/envs/topGO_env/bin/Rscript topGO_fixed.r
```

`topGO_env` is the only environment here with topGO installed. The batlow colours
are embedded as hex in the script, so `scico` is not required.

Set the contrast at the top of the script (`contrast <- "overall_clay_vs_sandy"`);
both directions are then processed in one run.

### Output — `overall_clay_vs_sandy/topgo_corrected/`

| file | contents |
|---|---|
| `GO_{up,down}_all_tested_terms.csv` | all 2,633 terms with counts, fold enrichment, classic p, BH-adjusted p, weight01 p |
| `GO_{up,down}_FDR0.05.csv` | terms passing FDR < 0.05 (67 up, 0 down) |
| `GO_{up,down}_weight01_p0.05.csv` | terms with weight01 p < 0.05 (28 up, 19 down) |
| `GO_up_classic_BH.{pdf,svg,png}` | main up-regulated figure |
| `GO_up_weight01.{pdf,svg,png}` | topology-corrected up-regulated figure |
| `GO_down_weight01.{pdf,svg,png}` | down-regulated, exploratory (not significant) |
| `GO_down_classic_rawp_NOT_significant.{pdf,svg,png}` | down-regulated, uncorrected p (not significant) |
