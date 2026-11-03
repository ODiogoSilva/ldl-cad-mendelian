#!/usr/bin/env Rscript
# Two-sample Mendelian randomization, both directions:
#   LDL -> CAD (forward) and CAD -> LDL (reverse control).
# Instruments: genome-wide significant (p < 5e-8), LD-clumped vs 1000G EUR.
# Estimators: IVW (multiplicative random effects), MR-Egger, weighted median.
# Diagnostics: Cochran/Rucker Q heterogeneity tests, Egger-intercept
# pleiotropy test, leave-one-out, mean instrument F statistic.

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", Ncpus = 4)
}
ensure_pkg("data.table")
ensure_pkg("jsonlite")
library(data.table)
library(jsonlite)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
ldl_f        <- args[1]
cad_f        <- args[2]
ldl_clump_f  <- args[3]
cad_clump_f  <- args[4]
outdir       <- args[5]
cad_n        <- as.numeric(args[6])
cad_ncase    <- as.numeric(args[7])
cad_ncontrol <- as.numeric(args[8])
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ldl <- fread(ldl_f)
cad <- fread(cad_f)
ldl_iv <- fread(ldl_clump_f)[["SNP"]]
cad_iv <- fread(cad_clump_f)[["SNP"]]
message("LDL instruments after clumping: ", length(ldl_iv))
message("CAD instruments after clumping: ", length(cad_iv))

comp <- c(A = "T", T = "A", C = "G", G = "C")

harmonise <- function(exp, out) {
  m <- merge(exp, out, by = "SNP", suffixes = c("_X", "_Y"))
  m <- m[!is.na(beta_X) & !is.na(beta_Y) & se_X > 0 & se_Y > 0]
  m[, `:=`(A1_Yc = unname(comp[A1_Y]), A2_Yc = unname(comp[A2_Y]))]
  ok_same  <- m$A1_X == m$A1_Y  & m$A2_X == m$A2_Y
  ok_swap  <- m$A1_X == m$A2_Y  & m$A2_X == m$A1_Y
  ok_samec <- m$A1_X == m$A1_Yc & m$A2_X == m$A2_Yc
  ok_swapc <- m$A1_X == m$A2_Yc & m$A2_X == m$A1_Yc
  keep <- ok_same | ok_swap | ok_samec | ok_swapc
  m <- m[keep]
  flipv <- (ok_swap | ok_swapc)[keep]
  m[flipv, `:=`(beta_Y = -beta_Y, EAF_Y = 1 - EAF_Y)]
  # Palindromic SNPs: keep only when allele frequencies are present and
  # clearly concordant (same minor allele in both studies).
  pal <- with(m, (A1_X == "A" & A2_X == "T") | (A1_X == "T" & A2_X == "A") |
                 (A1_X == "C" & A2_X == "G") | (A1_X == "G" & A2_X == "C"))
  concordant <- (!is.na(m$EAF_X) & !is.na(m$EAF_Y)) &
    ((m$EAF_X < 0.42 & m$EAF_Y < 0.42) | (m$EAF_X > 0.58 & m$EAF_Y > 0.58))
  m <- m[!(pal & !concordant)]
  m[, c("A1_Yc", "A2_Yc") := NULL]
  m
}

run_mr <- function(h, n_boot = 2000) {
  bX <- h$beta_X; bY <- h$beta_Y; sX <- h$se_X; sY <- h$se_Y
  n  <- length(bX)
  r  <- bY / bX
  sr <- sY / abs(bX)
  w  <- 1 / sr^2

  ## Inverse-variance weighted (multiplicative random effects)
  b_ivw  <- sum(w * r) / sum(w)
  se_fe  <- sqrt(1 / sum(w))
  Q      <- sum(w * (r - b_ivw)^2)
  df     <- n - 1
  phi    <- max(1, sqrt(Q / df))
  se_ivw <- se_fe * phi
  p_ivw  <- 2 * pnorm(abs(b_ivw / se_ivw), lower.tail = FALSE)
  Q_p    <- pchisq(Q, df = df, lower.tail = FALSE)

  ## MR-Egger (weighted regression with intercept; SEs inflated by residual sigma)
  wy  <- 1 / sY^2
  fit <- summary(lm(bY ~ bX, weights = wy))
  b_eg <- coef(fit)["bX", "Estimate"];          se_eg <- coef(fit)["bX", "Std. Error"]
  p_eg <- 2 * pt(abs(b_eg / se_eg), df = n - 2, lower.tail = FALSE)
  a_eg <- coef(fit)["(Intercept)", "Estimate"]; se_a  <- coef(fit)["(Intercept)", "Std. Error"]
  p_a  <- 2 * pt(abs(a_eg / se_a), df = n - 2, lower.tail = FALSE)
  Qe   <- sum(wy * (bY - (a_eg + b_eg * bX))^2)
  Qe_p <- pchisq(Qe, df = n - 2, lower.tail = FALSE)

  ## Weighted median (parametric bootstrap SE)
  wm_est <- function(rr, ss) {
    ww <- 1 / ss^2
    o  <- order(rr)
    cs <- cumsum(ww[o]) / sum(ww)
    rr[o][which(cs >= 0.5)[1]]
  }
  b_wm <- wm_est(r, sr)
  set.seed(42)
  boots <- replicate(n_boot, {
    i <- sample.int(n, n, replace = TRUE)
    wm_est(rnorm(n, mean = r[i], sd = sr[i]), sr[i])
  })
  se_wm <- sd(boots)
  p_wm  <- 2 * pnorm(abs(b_wm / se_wm), lower.tail = FALSE)

  ## Leave-one-out (IVW)
  loo <- vapply(seq_len(n), function(i) {
    wi <- w[-i]
    sum(wi * r[-i]) / sum(wi)
  }, numeric(1))

  list(
    n_instruments = n,
    mean_F        = mean((bX / sX)^2),
    estimates = list(
      IVW             = list(beta = b_ivw, se = se_ivw,
                             ci_lo = b_ivw - 1.96 * se_ivw, ci_hi = b_ivw + 1.96 * se_ivw, p = p_ivw),
      MR_Egger        = list(beta = b_eg,  se = se_eg,
                             ci_lo = b_eg - 1.96 * se_eg,  ci_hi = b_eg + 1.96 * se_eg,  p = p_eg),
      Weighted_median = list(beta = b_wm,  se = se_wm,
                             ci_lo = b_wm - 1.96 * se_wm,  ci_hi = b_wm + 1.96 * se_wm,  p = p_wm)
    ),
    heterogeneity = list(Q_ivw = Q, Q_df = df, Q_ivw_p = Q_p,
                         Q_egger = Qe, Q_egger_df = n - 2, Q_egger_p = Qe_p),
    pleiotropy = list(egger_intercept = a_eg, intercept_se = se_a, intercept_p = p_a),
    leave_one_out = list(min_beta = min(loo), max_beta = max(loo),
                         all_same_sign = (all(loo > 0) || all(loo < 0))),
    loo = data.table(SNP = h$SNP, loo_beta = loo)
  )
}

## ---------------- Forward: LDL -> CAD ----------------
exp1 <- ldl[SNP %in% ldl_iv]
h1   <- harmonise(exp1, cad)
message("Forward instruments after harmonization: ", nrow(h1))
fwrite(h1, file.path(outdir, "forward_harmonised.csv"))
res_fwd <- run_mr(h1)
fwrite(res_fwd$loo, file.path(outdir, "forward_leave_one_out.csv"))
loo_fwd <- res_fwd$loo
res_fwd$loo <- NULL
for (nm in names(res_fwd$estimates)) {
  e <- res_fwd$estimates[[nm]]
  res_fwd$estimates[[nm]]$OR       <- exp(e$beta)
  res_fwd$estimates[[nm]]$OR_ci_lo <- exp(e$ci_lo)
  res_fwd$estimates[[nm]]$OR_ci_hi <- exp(e$ci_hi)
}

## ---------------- Reverse: CAD -> LDL ----------------
exp2 <- cad[SNP %in% cad_iv]
h2   <- harmonise(exp2, ldl)
message("Reverse instruments after harmonization: ", nrow(h2))
fwrite(h2, file.path(outdir, "reverse_harmonised.csv"))
res_rev <- run_mr(h2)
fwrite(res_rev$loo, file.path(outdir, "reverse_leave_one_out.csv"))
res_rev$loo <- NULL

## ---------------- Plots ----------------
mk_plots <- function(h, res, tag, xlab, ylab) {
  d <- data.table(bX = h$beta_X, bY = h$beta_Y, sX = h$se_X, sY = h$se_Y)
  p1 <- ggplot(d, aes(x = bX, y = bY)) +
    geom_errorbarh(aes(xmin = bX - 1.96 * sX, xmax = bX + 1.96 * sX), alpha = 0.3) +
    geom_errorbar(aes(ymin = bY - 1.96 * sY, ymax = bY + 1.96 * sY), alpha = 0.3) +
    geom_point(size = 2) +
    geom_abline(slope = res$estimates$IVW$beta, intercept = 0, color = "#2166ac", linewidth = 1) +
    geom_abline(slope = res$estimates$MR_Egger$beta, intercept = res$pleiotropy$egger_intercept,
                color = "#b2182b", linewidth = 1, linetype = "dashed") +
    geom_abline(slope = res$estimates$Weighted_median$beta, intercept = 0,
                color = "#1b7837", linewidth = 1, linetype = "dotted") +
    labs(x = xlab, y = ylab, title = paste0("MR scatter: ", tag),
         caption = "Blue solid: IVW | red dashed: MR-Egger | green dotted: weighted median") +
    theme_bw()
  ggsave(file.path(outdir, paste0("scatter_", tag, ".png")), p1, width = 7, height = 6, dpi = 150)

  r  <- d$bY / d$bX
  sr <- d$sY / abs(d$bX)
  p2 <- ggplot(data.table(r = r, inv_se = 1 / sr), aes(x = r, y = inv_se)) +
    geom_point() +
    geom_vline(xintercept = res$estimates$IVW$beta, color = "#2166ac") +
    labs(x = "Wald ratio", y = "1 / SE(Wald ratio)", title = paste0("Funnel plot: ", tag)) +
    theme_bw()
  ggsave(file.path(outdir, paste0("funnel_", tag, ".png")), p2, width = 7, height = 6, dpi = 150)
}

mk_plots(h1, res_fwd, "forward", "LDL effect (beta, mmol/L per allele)", "CAD effect (log-OR per allele)")
mk_plots(h2, res_rev, "reverse", "CAD effect (log-OR per allele)", "LDL effect (beta, mmol/L per allele)")

# Leave-one-out forest (forward)
lo <- copy(loo_fwd)
lo[, SNP := factor(SNP, levels = SNP[order(loo_beta)])]
p3 <- ggplot(lo, aes(x = loo_beta, y = SNP)) +
  geom_point(size = 1) +
  geom_vline(xintercept = res_fwd$estimates$IVW$beta, color = "#b2182b") +
  labs(x = "IVW estimate leaving one SNP out (log-OR per mmol/L)", y = NULL,
       title = "Leave-one-out: LDL -> CAD") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 4))
ggsave(file.path(outdir, "loo_forward.png"), p3, width = 7, height = 10, dpi = 150, limitsize = FALSE)

## ---------------- JSON ----------------
out <- list(
  forward_LDL_to_CAD = res_fwd,
  reverse_CAD_to_LDL = res_rev,
  notes = list(
    cad_total_N = cad_n, cad_cases = cad_ncase, cad_controls = cad_ncontrol,
    forward_units = "log-OR(CAD) per 1 mmol/L increase in LDL cholesterol; OR columns give the odds ratio",
    reverse_units = "mmol/L change in LDL per 1 log-OR increase in CAD liability",
    reverse_expectation = "No (or negligible) causal effect of CAD liability on LDL cholesterol is expected; a clear effect here would indicate reverse causation or misspecification."
  )
)
write_json(out, file.path(outdir, "mr_results.json"), auto_unbox = TRUE,
           pretty = TRUE, digits = NA)
message("MR analysis done.")
