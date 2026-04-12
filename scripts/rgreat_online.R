#!/usr/bin/env Rscript

# Run online GREAT analysis for one BED file.
#
# To install dependencies, follow the guide on rGREAT repo.
# TODO: better idea to have the script install the package automatically,
# because whoever uses this might not know how to do it the better way.
# 
# The better way:
#   BiocManager::install(c("rGREAT", "GenomicRanges", "rtracklayer"), Ncpus = 128)
#   Yes, just use interact session and use a RM node, because compiling R packages
#   is painfully slow in RM-small and RM-shared. I learned this the hard way.
#
# Usage:
#   Rscript scripts/rgreat_online.R <bed_file> <genome> <output_dir>
#
# Example:
#   Rscript scripts/rgreat_online.R \
#     results/bedtools/cross_species_ep/human_open_in_mouse_promoters.bed \
#     mm10 \
#     results/rgreat/human_open_in_mouse_promoters

# check for arguments first
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scripts/rgreat_online.R <bed_file> <genome> <output_dir>")
}

bed_file <- args[1]
genome <- args[2]      # expected: "mm10" or "hg38"
output_dir <- args[3]

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rGREAT)
})

# do basic input checks

if (!file.exists(bed_file)) {
  stop(paste("BED file does not exist:", bed_file))
}

if (!(genome %in% c("mm10", "hg38"))) {
  stop("Genome must be one of: mm10, hg38")
}

# create the output directory for this one GREAT run.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# use the BED filename stem as the run name
# example:
#   human_open_in_mouse_promoters.bed
# becomes:
#   human_open_in_mouse_promoters
run_name <- tools::file_path_sans_ext(basename(bed_file))

# 
# read BED file
# 
# GREAT only needs genomic intervals, so we read only the first 3 columns:
#   chr, start, end
bed_df <- read.table(
  bed_file,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE
)

# panic if bed has less than 3 columns.
if (ncol(bed_df) < 3) {
  stop("BED file must have at least 3 columns: chr, start, end")
}

bed_df <- bed_df[, 1:3]
colnames(bed_df) <- c("chr", "start", "end")

# make sure start/end are numeric.
bed_df$start <- as.integer(bed_df$start)
bed_df$end <- as.integer(bed_df$end)

# panic if values are not numeric.
if (any(is.na(bed_df$start)) || any(is.na(bed_df$end))) {
  stop("BED start/end columns contain non-numeric values")
}

# sanity check but might not be neccessary: remove invalid rows where end <= start.
invalid_rows <- bed_df$end <= bed_df$start
if (any(invalid_rows)) {
  warning(paste("Removing", sum(invalid_rows), "invalid BED rows with end <= start"))
  bed_df <- bed_df[!invalid_rows, , drop = FALSE]
}

if (nrow(bed_df) == 0) {
  stop("No valid BED intervals remain after filtering")
}

# 
# convert BED to GRanges
# 
# BED uses 0-based, half-open coordinates.
# GRanges uses 1-based, closed coordinates.
# https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
gr <- GRanges(
  seqnames = bed_df$chr,
  ranges = IRanges(start = bed_df$start + 1, end = bed_df$end)
)

# Submit to an online job
# THIS IS WHORE THE MAGIC HAPPENS.
job <- submitGreatJob(
  gr,
  genome = genome
)

# retrieve all enrichment tables returned by GREAT.
tbls <- getEnrichmentTables(job)

# save the GREAT job object so it can be reloaded later in R.
saveRDS(job, file.path(output_dir, paste0(run_name, ".job.rds")))

# save each enrichment table as a CSV.
for (nm in names(tbls)) {
  safe_name <- gsub("[^A-Za-z0-9]+", "_", nm)
  out_csv <- file.path(output_dir, paste0(run_name, ".", safe_name, ".csv"))
  write.csv(tbls[[nm]], out_csv, row.names = FALSE)
}

# save a tiny metadata file for bookkeeping.
metadata_path <- file.path(output_dir, paste0(run_name, ".metadata.txt"))
writeLines(
  c(
    paste("bed_file:", bed_file),
    paste("genome:", genome),
    paste("run_name:", run_name),
    paste("n_regions:", length(gr))
  ),
  con = metadata_path
)

message("Finished GREAT online analysis for: ", bed_file)
message("Results written to: ", output_dir)

