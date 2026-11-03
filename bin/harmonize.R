#!/usr/bin/env Rscript
# Harmonize raw LDL (GLGC 2013) and CAD (CARDIoGRAMplusC4D 2015) summary
# statistics into one standard schema. Both files are already on the
# GRCh37/hg19 build (LDL via its SNP_hg19 column, CAD via bp_hg19); this step
# parses coordinates, applies QC and documents the builds explicitly.

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", Ncpus = 4)
}
ensure_pkg("data.table")
ensure_pkg("jsonlite")
library(data.table)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
ldl_file <- args[1]
cad_file <- args[2]
outdir   <- args[3]

valid_bases <- c("A", "C", "G", "T")

## ---------------- LDL exposure: GLGC / Willer et al. 2013 ----------------
message("Reading LDL file: ", ldl_file)
ldl <- fread(cmd = paste("zcat", shQuote(ldl_file)))
stopifnot(all(c("rsid", "A1", "A2", "beta", "se", "N", "P-value", "SNP_hg19") %in% names(ldl)))
setnames(ldl, c("rsid", "P-value", "Freq.A1.1000G.EUR"), c("SNP", "P", "EAF"))
parts <- tstrsplit(ldl[["SNP_hg19"]], ":", fixed = TRUE)
ldl[, `:=`(chr = gsub("chr", "", parts[[1]]),
           pos = suppressWarnings(as.integer(parts[[2]])))]
ldl <- ldl[, .(SNP,
               A1 = toupper(A1), A2 = toupper(A2),
               beta = as.numeric(beta), se = as.numeric(se), P = as.numeric(P),
               EAF = as.numeric(EAF), N = as.numeric(N), chr, pos)]
ldl <- ldl[complete.cases(ldl)]
ldl <- ldl[A1 %in% valid_bases & A2 %in% valid_bases]
ldl <- ldl[!duplicated(SNP)]
ldl <- ldl[P > 0 & se > 0]
fwrite(ldl, file.path(outdir, "ldl_std.tsv.gz"), sep = "\t")
fwrite(ldl[, .(SNP, P)], file.path(outdir, "ldl_snps_p.tsv"), sep = "\t")
message("LDL variants kept: ", nrow(ldl))

## ---------------- CAD outcome: CARDIoGRAMplusC4D 2015 ----------------
message("Reading CAD file: ", cad_file)
cad <- fread(cad_file)
stopifnot(all(c("markername", "chr", "bp_hg19", "effect_allele", "noneffect_allele",
                "beta", "se_dgc", "p_dgc") %in% names(cad)))
setnames(cad,
         c("markername", "chr", "bp_hg19", "effect_allele", "noneffect_allele",
           "effect_allele_freq", "beta", "se_dgc", "p_dgc"),
         c("SNP", "CHR", "BP", "A1_raw", "A2_raw", "EAF", "beta", "se", "P"))
cad <- cad[, .(SNP,
               A1 = toupper(A1_raw), A2 = toupper(A2_raw),
               beta = as.numeric(beta), se = as.numeric(se), P = as.numeric(P),
               EAF = suppressWarnings(as.numeric(EAF)), N = 184305,
               chr = as.character(CHR), pos = as.integer(BP))]
cad <- cad[complete.cases(cad[, .(SNP, A1, A2, beta, se, P, chr, pos)])]
cad <- cad[A1 %in% valid_bases & A2 %in% valid_bases]
cad <- cad[!duplicated(SNP)]
cad <- cad[P > 0 & se > 0]
fwrite(cad, file.path(outdir, "cad_std.tsv.gz"), sep = "\t")
fwrite(cad[, .(SNP, P)], file.path(outdir, "cad_snps_p.tsv"), sep = "\t")
message("CAD variants kept: ", nrow(cad))

info <- list(
  exposure = list(
    trait = "LDL cholesterol",
    study = "Global Lipids Genetics Consortium, Willer et al. 2013 (Nat Genet 45:1274-1283)",
    url   = "https://csg.sph.umich.edu/willer/public/lipids2013/jointGwasMc_LDL.txt.gz",
    sample_size_max = 188577,
    ancestry = "European",
    genome_build = "GRCh37/hg19 (from SNP_hg19 column)",
    variants = nrow(ldl)
  ),
  outcome = list(
    trait = "Coronary artery disease",
    study = "CARDIoGRAMplusC4D, Nikpay et al. 2015 (Nat Genet 47:1121-1130); GWAS Catalog GCST003116",
    url   = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt",
    cases = 60801, controls = 123504, sample_size = 184305,
    ancestry = "European",
    genome_build = "GRCh37/hg19 (bp_hg19 column)",
    variants = nrow(cad)
  ),
  ld_reference = list(
    panel = "1000 Genomes Phase 3, European (EUR) subset, PLINK format (MAGMA g1000_eur)",
    source = "Zenodo record 14142263 (g1000_eur.zip)",
    samples = 503, genome_build = "GRCh37/b37"
  )
)
write_json(info, file.path(outdir, "harmonise_info.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = NA)
message("Harmonization done.")
