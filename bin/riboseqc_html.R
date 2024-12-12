#!/usr/bin/env Rscript

# Load libraries ----------------------------------------------------------
message("Loading required libraries ...")
suppressPackageStartupMessages({
  library(RiboseQC)
  library(rmarkdown)
})

# Get variables from input ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_files <- args[1]
pandoc_dir <- args[2]
render_file <- args[3]
orfquant_prefix <- args[4]

# Define variables --------------------------------------------------------
rmarkdown::find_pandoc(dir = pandoc_dir)

input_files <- list.files(input_dir,
                          recursive = TRUE,
                          pattern = "_results_RiboseQC_all",
                          full.names = TRUE)
input_sample_names <- dirname(basename(input_files))

# Run script --------------------------------------------------------------

# get input and output file paths
input_files <- paste(normalizePath(dirname(input_files)),
                     basename(input_files), sep = "/")
output_dir <- paste(normalizePath(dirname(output_dir)),
                    basename(output_dir), sep = "/")

# set tmp folder outside R package
tmp_dir <- dirname(output_dir)

# create folder for rds objects and pdf figures
dir.create(paste0(output_dir, "/rds/"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(paste0(output_dir, "/pdf/"),
           recursive = TRUE, showWarnings = FALSE)
sink(file = paste0(output_dir, orfquant_prefix, "_report_text_output.txt"))

# render RMarkdown file > html report
knitclean <- knitr::knit_meta(class = NULL, clean = TRUE)
suppressWarnings(render(render_file,
                        params = list(input_files = input_files,
                                      input_sample_names = input_sample_names,
                                      output_fig_path = output_fig_path),
                        output_file = output_dir,
                        intermediates_dir = tmp_dir))
gici <- gc()
sink()
