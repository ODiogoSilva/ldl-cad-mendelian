#!/usr/bin/env Rscript
# Fine-map the PCSK9 locus for LDL cholesterol with SuSiE, using summary
# statistics (z-scores) and an in-sample LD matrix computed from the same
# 1000 Genomes EUR reference panel used for clumping.

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", Ncpus = 4)
}
ensure_pkg("data.table")
ensure_pkg("jsonlite")
ensure_pkg("BEDMatrix")
if (!requireNamespace("susieR", quietly = TRUE)) {
  ok <- FALSE
  if (requireNamespace("remotes", quietly = TRUE))
    ok <- tryCatch({ remotes::install_version("susieR", version = "0.12.35",
                                              repos = "https://cloud.r-project.org",
                                              upgrade = "never", Ncpus = 4); TRUE },
                   error = function(e) FALSE)
  if (!ok) install.packages("susieR", repos = "https://cloud.r-project.org", Ncpus = 4)
}
library(data.table)
library(jsonlite)
library(BEDMatrix)
library(susieR)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
ldl_f        <- args[1]
panel_prefix <- args[2]   # e.g. panel/g1000_eur (no extension)
outdir       <- args[3]
rchrom       <- args[4]
rstart       <- as.integer(args[5])
rend         <- as.integer(args[6])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ldl <- fread(ldl_f)
reg <- ldl[chr == rchrom & pos >= rstart & pos <= rend]
message("Variants in PCSK9 window: ", nrow(reg))

bim <- fread(paste0(panel_prefix, ".bim"),
             col.names = c("chr", "id", "cm", "pos_p", "A1_p", "A2_p"))
bim[, bim_idx := .I]
fam <- fread(paste0(panel_prefix, ".fam"),
             col.names = c("fid", "iid", "pid", "mid", "sex", "phe"))
message("Panel: ", nrow(bim), " variants, ", nrow(fam), " samples")

m <- merge(reg, bim[, .(id, bim_idx, pos_p, A1_p, A2_p)], by.x = "SNP", by.y = "id")
m <- m[(A1 == A1_p & A2 == A2_p) | (A1 == A2_p & A2 == A1_p)]
m[, flip := (A1 == A2_p)]
m <- m[!duplicated(SNP)]
message("Variants matched to panel with consistent alleles: ", nrow(m))

G <- BEDMatrix(panel_prefix)
Gm <- as.matrix(G[, m$bim_idx])
storage.mode(Gm) <- "numeric"

## Auto-detect allele dosage coding using the LDL file's 1000G-EUR allele
## frequencies (both refer to the same reference panel).
freq_obs <- colMeans(Gm, na.rm = TRUE) / 2
notflip  <- which(!m$flip)
probe    <- notflip[seq_len(min(300, length(notflip)))]
c_test   <- cor(freq_obs[probe], m$EAF[probe], use = "pairwise.complete.obs")
counts_A1 <- isTRUE(c_test > 0)
message("Allele coding check: cor(freq_obs, LDL EAF) = ", round(c_test, 3),
        " -> dosages count bim ", ifelse(counts_A1, "A1", "A2"))

## Recode to LDL effect-allele (A1) dosage
D <- Gm
for (j in seq_along(m$flip)) {
  if (counts_A1) {
    if (m$flip[j]) D[, j] <- 2 - Gm[, j]
  } else {
    if (!m$flip[j]) D[, j] <- 2 - Gm[, j]
  }
}
freq_eff <- colMeans(D, na.rm = TRUE) / 2
af_cor <- cor(freq_eff, m$EAF, use = "pairwise.complete.obs")
message("Effect-allele frequency concordance with LDL EAF: ", round(af_cor, 4))

## Impute missing genotypes with column means; drop zero-variance columns
for (j in seq_len(ncol(D))) {
  v <- D[, j]
  if (anyNA(v)) D[is.na(v), j] <- mean(v, na.rm = TRUE)
}
sds <- apply(D, 2, sd)
keepv <- sds > 0
D <- D[, keepv, drop = FALSE]
m <- m[keepv]
message("Variants after zero-variance filter: ", nrow(m))

z <- m$beta / m$se
R <- cor(D)
n_eff <- round(median(m$N, na.rm = TRUE))
message("Median sample size in region: ", n_eff)

fit <- tryCatch(
  susie_rss(z = z, R = R, n = n_eff, L = 10,
            estimate_residual_variance = TRUE, check_z = FALSE),
  error = function(e) {
    message("susie_rss with estimated residual variance failed: ", conditionMessage(e),
            " -> retrying with fixed residual variance")
    susie_rss(z = z, R = R, n = n_eff, L = 10,
              estimate_residual_variance = FALSE, check_z = FALSE)
  }
)

pip <- fit$pip
cs  <- fit$sets$cs
covg <- fit$sets$coverage
lead_i <- which.max(pip)

cs_list <- list()
if (!is.null(cs) && length(cs) > 0) {
  for (nm in names(cs)) {
    idx <- cs[[nm]]
    cs_list[[nm]] <- list(size = length(idx),
                          coverage = unname(covg[nm]),
                          variants = m$SNP[idx],
                          lead = m$SNP[idx[which.max(pip[idx])]],
                          max_pip = max(pip[idx]))
  }
}

res <- data.table(SNP = m$SNP, chr = m$chr, pos = m$pos, A1 = m$A1, A2 = m$A2,
                  EAF = m$EAF, beta = m$beta, se = m$se, P = m$P, N = m$N,
                  z = z, pip = pip)
res[, cs := ""]
if (length(cs_list) > 0)
  for (nm in names(cs_list)) res[SNP %in% cs_list[[nm]]$variants, cs := nm]
fwrite(res, file.path(outdir, "pcsk9_finemap.csv"))

info <- list(
  locus = list(name = "PCSK9", chr = rchrom, start = rstart, end = rend,
               genome_build = "GRCh37/hg19"),
  n_variants_input = nrow(reg),
  n_variants_analyzed = nrow(m),
  effect_allele_freq_concordance = af_cor,
  median_N = n_eff,
  converged = isTRUE(fit$converged),
  lead_variant = m$SNP[lead_i],
  lead_variant_pos = m$pos[lead_i],
  lead_variant_pip = pip[lead_i],
  lead_variant_P = m$P[lead_i],
  n_credible_sets = length(cs_list),
  credible_sets = cs_list
)
write_json(info, file.path(outdir, "finemap.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = NA)

## Regional association + PIP plot
res[, cs_label := ifelse(cs == "", "not in 95% CS", paste0("in ", cs))]
p <- ggplot(res, aes(x = pos / 1e6, y = -log10(P))) +
  geom_point(aes(color = pip, size = cs_label != "not in 95% CS")) +
  scale_color_gradient(low = "#9ecae1", high = "#cb181d", limits = c(0, 1)) +
  scale_size_manual(values = c("TRUE" = 3, "FALSE" = 1.2), guide = "none") +
  labs(x = paste0("Chromosome ", rchrom, " position (Mb, GRCh37)"),
       y = expression(-log[10](P)),
       color = "SuSiE PIP",
       title = "PCSK9 locus: LDL association and SuSiE posterior inclusion probabilities",
       caption = "Large points: variants in a 95% credible set") +
  theme_bw()
ggsave(file.path(outdir, "pcsk9_region.png"), p, width = 8, height = 6, dpi = 150)

message("Fine-mapping done. Lead: ", m$SNP[lead_i], " (PIP=", round(pip[lead_i], 3), ")")
