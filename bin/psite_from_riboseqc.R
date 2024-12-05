#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

riboseqc_file <- args[1]
package_install_loc <- args[2]

# Set location of installed packages for this R version
paths <- c(package_install_loc,.libPaths())
.libPaths(paths)

suppressPackageStartupMessages({
    library(RiboseQC)
})

# Get basename of the file
fname <- gsub(pattern = "_for_ORFquant", x = basename(riboseqc_file), replacement = "")

# ORFquant files are R object files
load(riboseqc_file)

# Extract p-sites from list of objects
p_sites <- data.frame(for_ORFquant$p_sites_uniq)

# Extract columns to create BED file
bed <- data.frame(p_sites$seqnames, p_sites$start, p_sites$end, ".", p_sites$score, p_sites$strand)
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

data.table::fwrite(bed, paste0(fname, "_psites.bed"), quote = F, sep = "\t", col.names = F, row.names = F)
