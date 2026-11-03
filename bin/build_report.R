#!/usr/bin/env Rscript
# Render the PDF + interactive HTML report and write the machine-readable
# JSON summary (also echoed to stdout for the job log).

ensure_pkg <- function(p) {
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org", Ncpus = 4)
}
ensure_pkg("rmarkdown")
ensure_pkg("jsonlite")
ensure_pkg("plotly")
library(rmarkdown)
library(jsonlite)

args <- commandArgs(trailingOnly = TRUE)
rmd    <- args[1]
mrdir  <- args[2]
fmdir  <- args[3]
outdir <- args[4]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

render(rmd, output_format = "pdf_document",
       output_file = "LDL_CAD_MR_report.pdf", output_dir = outdir,
       params = list(mrdir = mrdir, fmdir = fmdir))
render(rmd, output_format = "html_document",
       output_file = "LDL_CAD_MR_report.html", output_dir = outdir,
       params = list(mrdir = mrdir, fmdir = fmdir))

mr <- read_json(file.path(mrdir, "mr_results.json"), simplifyVector = FALSE)
fm <- read_json(file.path(fmdir, "finemap.json"), simplifyVector = FALSE)

summary <- list(
  pipeline = "ldl-cad-mendelian v1.0.0",
  analysis = "Two-sample Mendelian randomization: LDL -> CAD (forward) and CAD -> LDL (reverse control), with SuSiE fine-mapping of the PCSK9 locus",
  forward_LDL_to_CAD = mr$forward_LDL_to_CAD,
  reverse_CAD_to_LDL = mr$reverse_CAD_to_LDL,
  reverse_expectation = mr$notes$reverse_expectation,
  pcsk9_finemap = fm
)
js <- toJSON(summary, auto_unbox = TRUE, pretty = TRUE, digits = NA)
writeLines(js, file.path(outdir, "summary.json"))

cat(js)
cat("\n")
message("Report built: ", file.path(outdir, "LDL_CAD_MR_report.pdf"),
        " and .html, summary.json written.")
