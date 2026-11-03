# ldl-cad-mendelian

**Does LDL cholesterol causally raise the risk of coronary artery disease?**
A fully reproducible two-sample Mendelian randomization (MR) pipeline built with
Nextflow DSL2, using only public GWAS summary statistics (no individual-level
data). It runs the forward analysis (LDL -> CAD), the reverse direction as a
negative control (CAD -> LDL), SuSiE fine-mapping of the PCSK9 locus, and
produces a PDF report, an interactive HTML report and a machine-readable JSON
summary.

---

## 1. Analysis design

| Step | What happens | Process |
|---|---|---|
| Download | Fetches the LDL and CAD summary statistics and the 1000 Genomes EUR LD reference panel (with retries) | `DOWNLOAD` |
| Harmonize | Standardizes both GWAS files to one schema (rsID, alleles, beta, se, P, EAF, N, chr, pos), QC, documents genome builds | `HARMONIZE` |
| Clump | Selects independent genome-wide-significant instruments (p < 5e-8) with PLINK `--clump` vs 1000G EUR (r2 < 0.001, 1000 kb) | `CLUMP` |
| MR | Allele harmonization, then IVW / MR-Egger / weighted-median in **both directions**, plus heterogeneity (Cochran Q, Rucker Q'), pleiotropy (Egger intercept), leave-one-out and instrument F statistics | `MR` |
| Fine-map | SuSiE fine-mapping of PCSK9 (chr1:55.2-55.8 Mb, GRCh37) with in-sample LD from the same 1000G EUR panel | `FINEMAP` |
| Report | PDF report, interactive HTML report (plotly), and `summary.json` printed to the job log | `REPORT` |

## 2. Input data (all public, no individual-level data)

| Dataset | Role | Study / citation | Sample size | Ancestry | Genome build | URL |
|---|---|---|---|---|---|---|
| LDL cholesterol | Exposure | Global Lipids Genetics Consortium (GLGC), Willer et al. 2013, *Nat Genet* 45:1274-1283 | up to 188,577 | European | GRCh37/hg19 | https://csg.sph.umich.edu/willer/public/lipids2013/jointGwasMc_LDL.txt.gz |
| Coronary artery disease | Outcome | CARDIoGRAMplusC4D, Nikpay et al. 2015, *Nat Genet* 47:1121-1130; GWAS Catalog accession **GCST003116** | 60,801 cases / 123,504 controls (184,305) | European | GRCh37/hg19 | https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt |
| 1000 Genomes Phase 3, EUR | LD reference panel | MAGMA `g1000_eur` PLINK panel, mirrored at Zenodo record 14142263 | 503 | European | GRCh37/b37 | https://zenodo.org/api/records/14142263/files/g1000_eur.zip/content |

**Why these studies?** Both GWAS predate UK Biobank and are European-ancestry,
matching the EUR LD reference. The CAD study (GWAS Catalog GCST003116) and the
GLGC lipids study share some CHARGE cohorts (e.g. ARIC, CHS, FHS, MESA), so an
estimated **~20,000-30,000 participants (~10-15%) may overlap**. In two-sample
MR, overlap at this level with strong instruments (every instrument has F > 29;
mean F is reported) biases estimates only slightly, and toward the observational
association rather than manufacturing one.

## 3. Tools, versions and containers

Every process runs in a public, version-pinned container.

| Process | Container | Key tool(s) |
|---|---|---|
| `DOWNLOAD` | `alpine:3.20.3` | curl 8.x, unzip (apk) |
| `HARMONIZE` | `rocker/verse:4.3.2` | R 4.3.2; data.table, jsonlite |
| `CLUMP` | `quay.io/biocontainers/plink:1.90b6.21--h7b50bb2_7` | PLINK 1.90b6.21 |
| `MR` | `rocker/verse:4.3.2` | R 4.3.2; data.table, jsonlite, ggplot2 |
| `FINEMAP` | `rocker/verse:4.3.2` | R 4.3.2; data.table, BEDMatrix, susieR 0.12.35 |
| `REPORT` | `rocker/verse:4.3.2` | R 4.3.2; rmarkdown, knitr, plotly, TeX Live |

MR estimators (IVW multiplicative random effects, MR-Egger, weighted median),
heterogeneity, pleiotropy and leave-one-out are implemented in base R in
`bin/mr_analysis.R` (formulas as in Burgess et al. and Bowden et al.), so no
GitHub-only packages are required at runtime; `susieR` is the only non-trivial
CRAN install and is version-pinned.

## 4. Repository layout

```
main.nf               Nextflow DSL2 workflow (one process per step)
nextflow.config       Parameters, per-process containers/resources
bin/harmonize.R       Standardize both GWAS files; document builds
bin/mr_analysis.R     Harmonize alleles; MR both directions + diagnostics + plots
bin/finemap_pcsk9.R   SuSiE fine-mapping at PCSK9 with 1000G EUR LD
bin/build_report.R    Render PDF/HTML; write + echo summary.json
bin/report.Rmd        Report template (introduction/methods/results/conclusions)
```

## 5. How to rerun

### Locally (needs Nextflow >= 22.10 and Docker or Apptainer/Singularity)

```bash
git clone https://github.com/ODiogoSilva/ldl-cad-mendelian.git
cd ldl-cad-mendelian
nextflow run main.nf --outdir results
```

Parameters (all in `nextflow.config`, overridable with `--name value`):

| Parameter | Default | Meaning |
|---|---|---|
| `outdir` | `results` | Output directory |
| `ldl_url` / `cad_url` / `panel_url` | see section 2 | Input URLs |
| `pval_threshold` | `5e-8` | Instrument p-value threshold |
| `clump_kb` / `clump_r2` | `1000` / `0.001` | LD clumping window / r2 |
| `pcsk9_chr` / `pcsk9_start` / `pcsk9_end` | `1` / `55205470` / `55830520` | Fine-mapping window (GRCh37) |

### On CloudOS

1. Push this repository to GitHub (branch `main`).
2. In CloudOS: **New pipeline > Import from GitHub**, repo
   `ODiogoSilva/ldl-cad-mendelian`, branch `main`, main file `main.nf`.
3. Launch as a batch job with the workspace default queue and the newest
   DSL2-capable Nextflow version. No extra parameters are required
   (`params.outdir` defaults to `results`).
4. Results appear under `results/`:
   - `00_input/` - raw downloads + unpacked LD panel
   - `01_harmonised/` - standardized GWAS files, `harmonise_info.json`
   - `02_clump/` - PLINK clump outputs (instrument lists)
   - `03_mr/` - harmonized instrument tables, `mr_results.json`, scatter/funnel/LOO PNGs
   - `04_finemap/` - `pcsk9_finemap.csv`, `finemap.json`, regional plot
   - `05_report/` - **`LDL_CAD_MR_report.pdf`**, **`LDL_CAD_MR_report.html`** (interactive), **`summary.json`**

The final process also prints `summary.json` (instrument counts, every MR
estimate with 95% CI and p-value, heterogeneity and pleiotropy tests, the PCSK9
lead variant and credible sets) to its standard output; it is echoed into the
Nextflow console between `=====RESULTS_JSON=====` markers so it can be read
straight from the job log.

## 6. Methods notes

- **Instruments**: genome-wide significant (p < 5e-8), LD-clumped vs 1000G EUR
  (r2 < 0.001, 1000 kb window).
- **Allele harmonization**: effect alleles aligned across studies; strand flips
  resolved from allele letters; palindromic SNPs (A/T, C/G) kept only when
  allele frequencies are present and clearly concordant in both studies.
- **Estimators**: IVW with multiplicative random effects; MR-Egger (slope with
  intercept as directional-pleiotropy test); weighted median with parametric
  bootstrap SE.
- **Sensitivity**: Cochran's Q (IVW) and Rucker's Q' (Egger) heterogeneity
  tests; leave-one-out IVW; reverse-direction MR as negative control.
- **Fine-mapping**: SuSiE on LDL z-scores with in-sample LD from the same
  1000G EUR panel; genotype/allele orientation auto-validated against the LDL
  file's 1000G-EUR allele frequencies.

## 7. Limitations

- MR's exclusion-restriction assumption can be violated by pleiotropic lipid
  loci; Egger and weighted-median estimates mitigate but do not eliminate this.
- ~10-15% estimated sample overlap (see section 2); strong instruments keep the
  expected bias small.
- European-ancestry data only; results may not transfer to other ancestries.
- Fine-mapping LD comes from 503 reference samples; credible sets would sharpen
  with a larger matched panel.

## 8. Data reuse statement

All inputs are public GWAS summary statistics and public reference panels; no
individual-level data are used. Please cite the original studies (Willer et al.
2013; Nikpay et al. 2015; 1000 Genomes Phase 3) when reusing results.
