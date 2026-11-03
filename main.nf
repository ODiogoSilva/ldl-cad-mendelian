nextflow.enable.dsl = 2

/*
 * Two-sample Mendelian randomization of LDL cholesterol on coronary artery
 * disease, plus the reverse direction and PCSK9 fine-mapping.
 * Steps: download -> harmonize -> clump -> MR (both directions) ->
 *        fine-mapping -> report. One process per step.
 */

process DOWNLOAD {
  publishDir "${params.outdir}/00_input", mode: 'copy'

  output:
    path 'ldl_raw.txt.gz', emit: ldl
    path 'cad_raw.txt',    emit: cad
    path 'panel',          emit: panel, type: 'dir'

  script:
  '''
  set -eu
  apk add --no-cache curl unzip ca-certificates > /dev/null 2>&1
  update-ca-certificates 2>/dev/null || true

  echo "[download] LDL exposure: GLGC / Willer et al. 2013"
  curl -fSL --retry 8 --retry-delay 20 --retry-all-errors -C - \
      -o ldl_raw.txt.gz "!{params.ldl_url}"

  echo "[download] CAD outcome: CARDIoGRAMplusC4D 2015 (GWAS Catalog GCST003116)"
  curl -fSL --retry 8 --retry-delay 20 --retry-all-errors -C - \
      -o cad_raw.txt "!{params.cad_url}"

  echo "[download] LD reference panel: 1000 Genomes EUR (Zenodo 14142263)"
  curl -fSL --retry 8 --retry-delay 20 --retry-all-errors -C - \
      -o g1000_eur.zip "!{params.panel_url}"
  mkdir -p panel
  unzip -o g1000_eur.zip -d panel

  ls -lh ldl_raw.txt.gz cad_raw.txt
  ls -lh panel
  '''
}

process HARMONIZE {
  publishDir "${params.outdir}/01_harmonised", mode: 'copy'

  input:
    path ldl
    path cad
    path rscript

  output:
    path 'ldl_std.tsv.gz',      emit: ldl_std
    path 'cad_std.tsv.gz',      emit: cad_std
    path 'ldl_snps_p.tsv',      emit: ldl_p
    path 'cad_snps_p.tsv',      emit: cad_p
    path 'harmonise_info.json', emit: info

  script:
  """
  Rscript ${rscript} ${ldl} ${cad} .
  """
}

process CLUMP {
  publishDir "${params.outdir}/02_clump", mode: 'copy'

  input:
    path panel
    path ldl_p
    path cad_p

  output:
    path 'ldl.clumped', emit: ldl_clumped
    path 'cad.clumped', emit: cad_clumped

  script:
  """
  plink --bfile ${panel}/g1000_eur \
    --clump ${ldl_p} --clump-field P --clump-snp-field SNP \
    --clump-p1 ${params.pval_threshold} --clump-p2 1 \
    --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb} \
    --allow-no-sex --memory 12000 --out ldl

  plink --bfile ${panel}/g1000_eur \
    --clump ${cad_p} --clump-field P --clump-snp-field SNP \
    --clump-p1 ${params.pval_threshold} --clump-p2 1 \
    --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb} \
    --allow-no-sex --memory 12000 --out cad
  """
}

process MR {
  publishDir "${params.outdir}/03_mr", mode: 'copy'

  input:
    path ldl_std
    path cad_std
    path ldl_clumped
    path cad_clumped
    path rscript

  output:
    path 'mr_out', type: 'dir', emit: mr_dir

  script:
  """
  Rscript ${rscript} ${ldl_std} ${cad_std} ${ldl_clumped} ${cad_clumped} mr_out \
      ${params.cad_n} ${params.cad_ncase} ${params.cad_ncontrol}
  """
}

process FINEMAP {
  publishDir "${params.outdir}/04_finemap", mode: 'copy'

  input:
    path ldl_std
    path panel
    path rscript

  output:
    path 'fm_out', type: 'dir', emit: fm_dir

  script:
  """
  Rscript ${rscript} ${ldl_std} ${panel}/g1000_eur fm_out \
      ${params.pcsk9_chr} ${params.pcsk9_start} ${params.pcsk9_end}
  """
}

process REPORT {
  publishDir "${params.outdir}/05_report", mode: 'copy'

  input:
    path mr_dir
    path fm_dir
    path rscript
    path rmd

  output:
    path 'report_out', type: 'dir', emit: report_dir
    path 'report_out/summary.json', emit: summary

  script:
  """
  Rscript ${rscript} ${rmd} ${mr_dir} ${fm_dir} report_out
  """
}

workflow {
  DOWNLOAD()
  HARMONIZE(DOWNLOAD.out.ldl, DOWNLOAD.out.cad, file("${projectDir}/bin/harmonize.R"))
  CLUMP(DOWNLOAD.out.panel, HARMONIZE.out.ldl_p, HARMONIZE.out.cad_p)
  MR(HARMONIZE.out.ldl_std, HARMONIZE.out.cad_std, CLUMP.out.ldl_clumped, CLUMP.out.cad_clumped,
     file("${projectDir}/bin/mr_analysis.R"))
  FINEMAP(HARMONIZE.out.ldl_std, DOWNLOAD.out.panel, file("${projectDir}/bin/finemap_pcsk9.R"))
  REPORT(MR.out.mr_dir, FINEMAP.out.fm_dir,
         file("${projectDir}/bin/build_report.R"), file("${projectDir}/bin/report.Rmd"))

  // Echo the machine-readable summary into the Nextflow console so it is
  // visible in the job log.
  REPORT.out.summary.view { f ->
    "\n=====RESULTS_JSON=====\n" + f.text + "\n=====END_RESULTS_JSON=====\n"
  }
}
