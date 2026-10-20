#!/usr/bin/env Rscript

# ============================================================
#  ShortReadsPhaser
#
#  Usage: Rscript run_pipeline.R --config config.yaml
#         Rscript run_pipeline.R --config config.yaml --redo
#         Rscript run_pipeline.R --config config.yaml --thread 8
# ============================================================

# --- Load libraries and function ----------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))

# create a function to retrieve parameter value with default
f_get_param <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

# create a function to check that a file or directory exists
f_check_exists <- function(path, label, is_dir=FALSE) {
  path  <- path.expand(path)

  is_ok <- if (is_dir) dir.exists(path) else file.exists(path)
  if (!is_ok) {
    stop(paste0(label, " not found: ", path))
  }
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="Path to YAML config file [required]", metavar="FILE"),
  make_option(c("-r", "--redo"), action="store_true", default=FALSE,
              help="Re-run all analyses and override previous results"),
  make_option(c("-t", "--thread"), type="integer", default=NULL,
              help="Number of threads (overrides the config file)", metavar="N")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if config file is invalid
if (is.null(opt$config)) {
  stop("-c/--config is required.\n
        Usage: Rscript run_pipeline.R --config config.yaml")
}

if (!file.exists(opt$config)) {
  stop(paste("Config file not found:", opt$config))
}

# --- Load and validate config --------------------------------
cfg <- yaml::read_yaml(opt$config)

# set required parameters
required_fields <- c("codedir", "outdir",
                     "dir_backbone_alignment", "dir_backbone_genetrees",
                     "file_backbone_genetree_models", "file_taxonomy_metadata", "file_backbone_speciestree",
                     "sample", "dir_fastq_gz", "extraction_tool", "phasing_tool")

# check if required fields are missing
missing <- setdiff(required_fields, names(cfg))
if (length(missing) > 0) {
  stop(paste("Missing required config fields:", paste(missing, collapse=", ")))
}

# check if required fields are empty
empty <- required_fields[sapply(required_fields, function(x) { is.null(cfg[[x]]) || identical(cfg[[x]], "") })]
if (length(empty) > 0) {
  stop(paste("Empty required config fields:", paste(empty, collapse=", ")))
}

# check if input files and directories exist
f_check_exists(cfg$dir_backbone_alignment, "dir_backbone_alignment", is_dir=TRUE)
f_check_exists(cfg$dir_backbone_genetrees, "dir_backbone_genetrees", is_dir=TRUE)
f_check_exists(cfg$file_backbone_genetree_models, "file_backbone_genetree_models")
f_check_exists(cfg$file_taxonomy_metadata, "file_taxonomy_metadata")
f_check_exists(cfg$file_backbone_speciestree, "file_backbone_speciestree")

# check if extraction_tool is valid
if (!tolower(cfg$extraction_tool) %in% c("captus", "hybpiper")) {
    stop(paste("Invalid extraction_tool:", cfg$extraction_tool, "(must be 'captus' or 'hybpiper')"))
}

# check if phasing_tool is valid
if (!tolower(cfg$phasing_tool) %in% c("whatshap", "hapcut2")) {
    stop(paste("Invalid phasing_tool:", cfg$phasing_tool, "(must be 'whatshap' or 'hapcut2')"))
}

# check if the input directories exist
if (tolower(cfg$extraction_tool) == "captus") {
  f_check_exists(cfg$file_query_matches_flanked, "file_query_matches_flanked")
  f_check_exists(cfg$file_query_matches_cds, "file_query_matches_cds")
} else if (tolower(cfg$extraction_tool) == "hybpiper") {
  f_check_exists(cfg$dir_hybpiper_query_output, "dir_hybpiper_query_output", is_dir=TRUE)
}

# check if the paired-end reads follow the expected naming (<sample>_R1/_R2.fastq.gz)
for (read in c("_R1.fastq.gz", "_R2.fastq.gz")) {
  f_check_exists(path  = file.path(cfg$dir_fastq_gz, paste0(cfg$sample, read)),
                 label = paste0(cfg$sample, read))
}

# check if prefix is set, otherwise use the sample name
if (is.null(cfg$prefix) || cfg$prefix == "") {
  cfg$prefix <- cfg$sample
}

# --- Apply CLI overrides -------------------------------------
if (opt$redo) {
  message("Note: --redo flag set via CLI, overriding config.")
  cfg$redo <- TRUE
}

if (!is.null(opt$thread)) {
  message("Note: --thread set via CLI, overriding config.")
  cfg$thread <- opt$thread
}

# --- Map config to rmarkdown params --------------------------
render_params <- list(
  codedir                       = cfg$codedir,
  prefix                        = cfg$prefix,
  outdir                        = cfg$outdir,
  thread                        = as.integer(f_get_param(cfg$thread, 1)),
  redo                          = as.logical(f_get_param(cfg$redo, FALSE)),

  exe_amas                      = f_get_param(cfg$exe_amas, "AMAS.py"),
  exe_bcftools                  = f_get_param(cfg$exe_bcftools, "bcftools"),
  exe_bwa                       = f_get_param(cfg$exe_bwa, "bwa"),
  exe_epa_ng                    = f_get_param(cfg$exe_epa_ng, "epa-ng"),
  exe_extracthairs              = f_get_param(cfg$exe_extracthairs, "extractHAIRS"),
  exe_gappa                     = f_get_param(cfg$exe_gappa, "gappa"),
  exe_hapcut2                   = f_get_param(cfg$exe_hapcut2, "HAPCUT2"),
  exe_mafft                     = f_get_param(cfg$exe_mafft, "mafft"),
  exe_iqtree                    = f_get_param(cfg$exe_iqtree, "iqtree2"),
  exe_samtools                  = f_get_param(cfg$exe_samtools, "samtools"),
  exe_whatshap                  = f_get_param(cfg$exe_whatshap, "whatshap"),

  dir_backbone_alignment        = cfg$dir_backbone_alignment,
  dir_backbone_genetrees        = cfg$dir_backbone_genetrees,
  file_backbone_genetree_models = cfg$file_backbone_genetree_models,
  file_taxonomy_metadata        = cfg$file_taxonomy_metadata,
  file_backbone_speciestree     = cfg$file_backbone_speciestree,

  sample                        = cfg$sample,
  dir_fastq_gz                  = cfg$dir_fastq_gz,
  extraction_tool               = tolower(cfg$extraction_tool),
  phasing_tool                  = tolower(cfg$phasing_tool),

  file_query_matches_flanked    = cfg$file_query_matches_flanked,
  file_query_matches_cds        = cfg$file_query_matches_cds,

  dir_hybpiper_query_output     = cfg$dir_hybpiper_query_output,

  min_depth                     = as.numeric(f_get_param(cfg$min_depth, 10)),
  max_depth                     = as.numeric(f_get_param(cfg$max_depth, 80)),
  min_coverage                  = as.numeric(f_get_param(cfg$min_coverage, 0.9)),
  min_identity                  = as.numeric(f_get_param(cfg$min_identity, 0.9)),
  taxonomy_rank                 = as.integer(f_get_param(cfg$taxonomy_rank, 1)),
  min_alwr                      = as.numeric(f_get_param(cfg$min_alwr, 0.8)),
  min_consensus_thresh          = as.numeric(f_get_param(cfg$min_consensus_thresh, 1.0))
)

# --- Run ShortReadsPhaser ------------------------------------
rmd_path <- file.path(path.expand(render_params$codedir), "codes", "1_main.Rmd")
if (!file.exists(rmd_path)) {
  stop(paste("1_main.Rmd not found:", rmd_path))
}

message("Starting ShortReadsPhaser pipeline...")
message("  Config:  ", opt$config)
message("  Sample:  ", render_params$sample)
message("  Prefix:  ", render_params$prefix)
message("  Output:  ", render_params$outdir)
message("  Threads: ", render_params$thread)

# render the Rmarkdown file
rmarkdown::render(
  input       = rmd_path,
  params      = render_params,
  output_file = paste0(render_params$prefix, "_report.html"),
  output_dir  = file.path(path.expand(render_params$outdir), render_params$prefix),
  quiet       = FALSE
)

message("Done. Report: ",
        file.path(path.expand(render_params$outdir), render_params$prefix,
                  paste0(render_params$prefix, "_report.html")))