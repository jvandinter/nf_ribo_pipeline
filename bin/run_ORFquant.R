#!/usr/bin/env Rscript

message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(ORFquant)
  library(RiboseQC)
  library(rmarkdown)
})

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
for_ORFquant_file <- args[1]
pool_id <- args[2]
rannot <- args[3]
cpu <- args[4]
pandoc_dir <- args[5]
annotation_package <- args[6]
package_install_loc <- args[7]
outdir <- args[8]

paths <- c(package_install_loc, .libPaths())
.libPaths(paths)

# Define variables --------------------------------------------------------
html_report <- paste0(pool_id,"_report.html")
ORFquant_output_file <- paste0(pool_id, "_final_ORFquant_results")
ORFquant_plot_data <- paste0(ORFquant_output_file, "_plots/",
                             pool_id, "_ORFquant_plots_RData")
find_pandoc(dir = pandoc_dir)

# Define functions --------------------------------------------------------
ORFquant_analysis <- function(for_ORFquant_file,
                              ORFquant_output_file,
                              ORFquant_plot_data,
                              rannot,
                              cpu,
                              outdir,
                              html_report,
                              sample_id) {

  print(for_ORFquant_file)

  if (!file.exists(ORFquant_output_file)) {

    run_ORFquant(for_ORFquant_file = for_ORFquant_file,
      annotation_file = rannot,
      n_cores = cpu,
      prefix = outdir,
      gene_name = NA,
      gene_id = NA,
      genomic_region = NA,
      write_temp_files = TRUE,
      write_GTF_file = TRUE,
      write_protein_fasta = TRUE,
      interactive = TRUE,
      stn.orf_find.all_starts = TRUE,
      stn.orf_find.nostarts = FALSE,
      stn.orf_find.start_sel_cutoff = NA,
      stn.orf_find.start_sel_cutoff_ave = 0.5,
      stn.orf_find.cutoff_fr_ave = 0.5,
      stn.orf_quant.cutoff_cums = NA,
      unique_reads_only = FALSE,
      canonical_start_only = TRUE,
      stn.orf_quant.scaling = "total_Psites"
    )
  } else {
    print("ORFquant file already exists")
  }

  plot_ORFquant_results(for_ORFquant_file = for_ORFquant_file,
    ORFquant_output_file = ORFquant_output_file,
    annotation_file = rannot
  )

  create_ORFquant_html_report(input_files = ORFquant_plot_data,
    input_sample_name = sample_id,
    output_file = html_report
  )
}

# Run script --------------------------------------------------------------
if (!require(basename(annotation_package), character.only = TRUE)) {
  install.packages(annotation_package,
                   repos = NULL, type = "source")
}

ORFquant_analysis(for_ORFquant_file = for_ORFquant_file,
                  ORFquant_output_file = ORFquant_output_file,
                  ORFquant_plot_data = ORFquant_plot_data,
                  rannot = rannot,
                  cpu = cpu,
                  html_report = html_report,
                  sample_id = sample_id)
