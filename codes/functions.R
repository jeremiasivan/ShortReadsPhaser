# functions for running ShortReadsPhaser

# function: extract all CAPTUS top hits
f_extract_captus_best_hits <- function(fn_captus_matches, fn_out) {
    # open the CAPTUS output
    seq <- Biostrings::readBStringSet(fn_captus_matches)

    # extract top hits for all loci
    best_hits <- seq[grepl("\\[hit=00\\]", names(seq))]

    # update sequence headers
    names(best_hits) <- sapply(names(best_hits), function(x) {
        locus <- unlist(strsplit(x, split="__"))[2]
        locus <- unlist(strsplit(locus, split=" "))[1]
        locus
    })

    # save output file
    Biostrings::writeXStringSet(best_hits, filepath=fn_out)
}

# function: create individual files for each locus
f_split_captus_best_hits <- function(fn_captus_best_hits, outdir) {
    # open the CAPTUS output
    seq <- Biostrings::readBStringSet(fn_captus_best_hits)

    # extract individual loci
    for (locus in names(seq)) {
        subseq <- seq[locus]

        # save the FASTA sequence
        fn_out <- file.path(outdir, paste0(locus, ".fna"))
        Biostrings::writeXStringSet(subseq, filepath=fn_out)
    }
}

# function: run BWA-MEM and index the BAM file
f_bwa_mem <-  function(fn_target_loci, fn_fastq_r1, fn_fastq_r2, fn_bam_sort, fn_bam_dedup, fn_depth, thread, exe_bwa, exe_samtools) {
    # index files
    system(paste(exe_samtools, "faidx", fn_target_loci))
    system(paste(exe_bwa, "index", fn_target_loci))

    # set Samtools thread
    samtools_thread <- paste0("-@", thread)

    # run BWA-MEM and sort the BAM file
    cmd_bwa <- paste(exe_bwa, "mem",
                     "-t", thread,
                     "-R '@RG\\tID:sample1\\tSM:sample1'",
                     fn_target_loci,
                     fn_fastq_r1, fn_fastq_r2,
                     "|", exe_samtools, "fixmate -m", samtools_thread, "- -",
                     "|", exe_samtools, "sort", samtools_thread, "-o", fn_bam_sort, "-")
    system(cmd_bwa)

    # deduplicate reads
    cmd_dedup <- paste(exe_samtools, "markdup", samtools_thread, "-r", fn_bam_sort, "-",
                       "|", exe_samtools, "view", "-b -q 20 -F 0x904", "-o", fn_bam_dedup, "-",
                       "&&", exe_samtools, "index", fn_bam_dedup)
    system(cmd_dedup)

    # calculate per-locus depth
    cmd_depth <- paste(exe_samtools, "depth", "-a", fn_bam_dedup, ">", fn_depth)
    system(cmd_depth)
}

# function: variant calling
f_variant_calling <- function(fn_target_loci, fn_bam, fn_vcf_gz, fn_vcf_gz_filtered, min_depth, max_depth, thread, exe_bcftools) {
    # do variant calling
    cmd_vcf <- paste(exe_bcftools, "mpileup",
                     "-f", fn_target_loci,
                     "-q 20 -Q 20 -a AD,DP,SP",
                     "--threads", thread,
                     fn_bam,
                     "|", exe_bcftools, "call", "-mv --ploidy 2 -Oz", "-o", fn_vcf_gz)
    system(cmd_vcf)
    system(paste(exe_bcftools, "index -t", fn_vcf_gz))

    # filter variants
    cmd_vcf_filter <- paste(exe_bcftools, "view",
                            "-m2 -M2 -i 'QUAL>=30 &&", paste0("INFO/DP>=", min_depth), "&&", paste0("INFO/DP<=", max_depth, "'"),
                            fn_vcf_gz,
                            "|", exe_bcftools, "filter", "-e 'GT=\"het\" && (FMT/AD[0:1] < 0.2*FMT/DP || FMT/AD[0:1] > 0.8*FMT/DP)'", "-s LOWAB",
                            "|", exe_bcftools, "view", "-f PASS,. -Oz", "-o", fn_vcf_gz_filtered)
    system(cmd_vcf_filter)
    system(paste(exe_bcftools, "index -t", fn_vcf_gz_filtered))
}

# function: phasing
f_whatshap <- function(fn_target_loci, fn_vcf_gz_filtered, fn_bam, fn_phased_vcf, fn_blocks_whatshap, exe_whatshap) {
    cmd_whatshap <- paste(exe_whatshap, "phase",
                          "--reference", fn_target_loci,
                          "--indels --distrust-genotypes --tag PS",
                          "-o", fn_phased_vcf,
                          fn_vcf_gz_filtered, fn_bam)
    system(cmd_whatshap)

    # compressed the output file
    system(paste("bgzip -f", fn_phased_vcf, "&&", "tabix -f -p vcf", paste0(fn_phased_vcf, ".gz")))

    # extract variant coordinates
    system(paste(exe_whatshap, "stats --block-list", fn_blocks_whatshap, paste0(fn_phased_vcf, ".gz")))
}

# function: generate haplotypes
f_generate_haplotypes <- function(fn_target_loci, fn_bed, fn_snps_vcf, fn_hap1, fn_hap2, exe_bcftools) {
    system(paste(exe_bcftools, "consensus -f", fn_target_loci, "-s sample1 -H 1pIu", "-m", fn_bed, fn_snps_vcf, ">", fn_hap1))
    system(paste(exe_bcftools, "consensus -f", fn_target_loci, "-s sample1 -H 2pIu", "-m", fn_bed, fn_snps_vcf, ">", fn_hap2))
}

# function: extract exon coordinates from a flanked sequence (source: Claude)
f_locate_exons <- function(fn_cds, fn_flanked, min_coverage, min_identity) {
    # open the sequences
    cds <- as.character(Biostrings::readBStringSet(fn_cds)[[1]])
    ref <- as.character(Biostrings::readBStringSet(fn_flanked)[[1]])

    # remove padded Ns
    cds <- toupper(gsub("n", "", cds, fixed=TRUE))
    if (nchar(cds) == 0) {
        return(NULL)
    }

    # extract the first and last positions across all exons in the flanked sequence
    exon_stretch <- gregexpr("[[:upper:]]+", ref)[[1]]
    if (exon_stretch[1] == -1) {
        return(NULL)
    }

    # extract the start and end positions of all exons
    exon_start <- as.integer(exon_stretch)
    exon_end   <- exon_start + attr(exon_stretch, "match.length") - 1L

    # concatenate all exons
    exon_concat <- toupper(paste(substring(ref, exon_start, exon_end), collapse=""))
    exon_coords <- unlist(Map(seq, exon_start, exon_end))
    
    # compare the CDS and concatenated exons
    if (exon_concat == cds) {
        positions <- exon_coords
        identity  <- 1
    } else {
        aln <- pwalign::pairwiseAlignment(Biostrings::BStringSet(cds), Biostrings::BStringSet(exon_concat), type="global")
        aln_cds   <- unlist(strsplit(as.character(pwalign::pattern(aln)), ""))
        aln_ref   <- unlist(strsplit(as.character(pwalign::subject(aln)), ""))
 
        # extract aligned positions and calculate identity
        exon_idx  <- cumsum(aln_ref != "-") + pwalign::start(pwalign::subject(aln)) - 1
        mapped    <- aln_cds != "-" & aln_ref != "-"
        positions <- exon_coords[exon_idx[mapped]]

        identity  <- mean(aln_cds[mapped] == aln_ref[mapped])
        if (length(positions) < min_coverage * nchar(cds) || identity < min_identity) {
            return(NULL)
        }
    }
 
    # collapse consecutive positions into intervals
    breaks    <- which(diff(positions) != 1)
    df_output <- data.table::data.table(strand="+", identity=round(identity, 4), unmapped=nchar(cds) - length(positions),
                                        start=positions[c(1, breaks + 1)], end=positions[c(breaks, length(positions))])
 
    return (df_output)
}

# function: extract exons out of a haplotype sequence (source: Claude)
f_extract_cds <- function(fn_hap, exons) {
    # open the haplotype sequence
    hap <- Biostrings::readBStringSet(fn_hap)
    seq <- toupper(as.character(hap[[1]]))
    if (exons$strand[1] == "-") {
        seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq)))
    }

    # extract the coding sequence
    cds <- Biostrings::BStringSet(setNames(paste(substring(seq, exons$start, exons$end), collapse=""), names(hap)[1]))
    return(cds)
}

# run MAFFT --add
f_mafft_add <- function(fn_ref, fn_sample, fn_out, exe_mafft) {
  cmd_mafft <- paste(exe_mafft, "--auto --add", fn_sample, "--keeplength --adjustdirection", fn_ref, ">", fn_out)
  system(cmd_mafft)
}

# substitute IQ-TREE2 models to EPA-NG models (source: Claude)
f_iqtree2epa_ng_model <- function(model) {
    map <- c(JC69="JC", K2P="K80", HKY85="HKY",
             TN="TN93", TrN="TN93", TNe="TN93ef",
             K3P="K81", TPM1="K81",
             K81u="K81uf", K3Pu="K81uf", TPM1u="K81uf",
             TPM2u="TPM2uf", TPM3u="TPM3uf",
             TIM="TIM1uf", TIMe="TIM1", TIM2="TIM2uf", TIM2e="TIM2", TIM3="TIM3uf", TIM3e="TIM3",
             TVMe="TVMef")

    # extract substitution model and modifiers
    matrix <- sub("\\+.*", "", model)
    modifiers <- sub("^[^+]*", "", model)

    # update the substitution model if it is in the map
    if (matrix %in% names(map)) {
        matrix <- map[[matrix]]
    }

    model <- paste0(matrix, modifiers)
    return(model)
}

# run EPA-NG
f_epa_ng <- function(fn_ref, fn_query, fn_tree, model, outdir, fn_log, exe_epa_ng) {
  cmd_epa_ng <- paste(exe_epa_ng,
                      "--ref-msa", fn_ref,
                      "--tree", fn_tree,
                      "--model", model,
                      "--query", fn_query,
                      "-w", outdir, "--redo", ">>", fn_log)
  system(cmd_epa_ng)
}

# run gappa examine assign
f_gappa_assign <- function(fn_jplace, fn_taxon, outdir, log_file, exe_gappa) {
    cmd_gappa <- paste(exe_gappa, "examine assign",
                       "--jplace-path", fn_jplace,
                       "--taxon-file", fn_taxon,
                       "--per-query-results --allow-file-overwriting",
                       "--out-dir", outdir,
                       ">>", log_file)
    system(cmd_gappa)
}

# function: read gappa assignments
f_read_gappa_assignments <- function(dir_epa_ng, tax_rank) {
    # list output files from gappa examine assign
    fn_per_query <- list.files(dir_epa_ng, pattern="^per_query\\.tsv$", recursive=TRUE, full.names=TRUE)

    # combine all assignments into a data.frame
    ls_assign <- lapply(fn_per_query, function(file) {
        if (file.size(file) == 0) {
            return(NULL)
        }

        # open the file
        df <- data.table::fread(file)
        df$locus <- basename(dirname(dirname(file)))
        df$aLWR  <- as.numeric(df$aLWR)
        df$depth <- lengths(strsplit(df$taxopath, ";"))

        # extract the rank of interest and remove duplicates
        df <- df %>%
                filter(depth == tax_rank) %>%
                group_by(name) %>%
                slice_max(aLWR, n=1, with_ties=FALSE) %>%
                ungroup()

        df$lineage <- sapply(strsplit(df$taxopath, ";"), function(path) { path[tax_rank] })
        df
    })
    df_assign <- data.table::rbindlist(ls_assign)

    # extract information from the sequence name
    df_assign$PS      <- gsub(".*_PS([^_]+)_h[12]$", "\\1", df_assign$name)
    df_assign$hap     <- gsub(".*_(h[12])$", "\\1", df_assign$name)

    # subset the columns of interest
    df_assign <- df_assign %>% select(c("locus", "PS", "hap", "lineage", "aLWR"))

    return(df_assign)
}

# function: find the two parental lineages for every block (source: Claude)
f_assign_blocks <- function(df_assign, min_alwr, input_lineage_A=NULL, input_lineage_B=NULL) {
    # extract unique blocks
    df_blocks <- df_assign %>%
                    select("locus", "PS") %>%
                    unique() %>%
                    mutate(h1_lineage="none", h2_lineage="none", h1_aLWR=0, h2_aLWR=0)

    # iterate over each block
    for (i in 1:nrow(df_blocks)) {
        # extract the rows corresponding to the current block
        r1 <- df_assign[df_assign$locus==df_blocks$locus[i] & df_assign$PS==df_blocks$PS[i] & df_assign$hap == "h1", ]
        r2 <- df_assign[df_assign$locus==df_blocks$locus[i] & df_assign$PS==df_blocks$PS[i] & df_assign$hap == "h2", ]

        # assign the lineages and aLWR values to the block
        if (nrow(r1) > 0) {
            df_blocks$h1_lineage[i] <- r1$lineage[1]
            df_blocks$h1_aLWR[i]    <- r1$aLWR[1]
        }

        if (nrow(r2) > 0) {
            df_blocks$h2_lineage[i] <- r2$lineage[1]
            df_blocks$h2_aLWR[i]    <- r2$aLWR[1]
        }
    }

    # check which blocks have confident assignments
    is_conf1 <- df_blocks$h1_aLWR >= min_alwr
    is_conf2 <- df_blocks$h2_aLWR >= min_alwr
    is_both  <- is_conf1 & is_conf2

    # pair the lineages of the two haplotypes for each block
    pair <- paste(pmin(df_blocks$h1_lineage, df_blocks$h2_lineage),
                  pmax(df_blocks$h1_lineage, df_blocks$h2_lineage), sep=" | ")

    # calculate the number of blocks for each pair of lineages
    df_pairs <- data.table::data.table(sort(table(pair[is_both]), decreasing=TRUE))
    colnames(df_pairs) <- c("pair", "n_blocks")

    # extract the two parental lineages from the most common pair
    most_common_pair <- unlist(strsplit(df_pairs$pair[1], split=" | ", fixed=TRUE))
    lineage_A <- ifelse(is.null(input_lineage_A), most_common_pair[1], input_lineage_A)
    lineage_B <- ifelse(is.null(input_lineage_B), most_common_pair[2], input_lineage_B)

    # assign decisions per block
    df_blocks$decision <- "ambiguous"
    df_blocks$tier     <- "none"

    for (i in 1:nrow(df_blocks)) {
        l1 <- df_blocks$h1_lineage[i]
        l2 <- df_blocks$h2_lineage[i]

        if (is_both[i] && l1 == lineage_A && l2 == lineage_B) {
            df_blocks$decision[i] <- "keep"
            df_blocks$tier[i]     <- "both"
        } else if (is_both[i] && l1 == lineage_B && l2 == lineage_A) {
            df_blocks$decision[i] <- "flip"
            df_blocks$tier[i]     <- "both"
        } else if (is_both[i] && l1 == l2) {
            df_blocks$decision[i] <- "same"
        } else if (is_conf1[i] && l1 %in% c(lineage_A, lineage_B) && !(l2 %in% c(lineage_A, lineage_B))) {
            df_blocks$decision[i] <- ifelse(l1 == lineage_A, "keep", "flip")
            df_blocks$tier[i]     <- "single"
        } else if (is_conf2[i] && l2 %in% c(lineage_A, lineage_B) && !(l1 %in% c(lineage_A, lineage_B))) {
            df_blocks$decision[i] <- ifelse(l2 == lineage_A, "flip", "keep")
            df_blocks$tier[i]     <- "single"
        }
    }

    return(list(blocks=df_blocks, pairs=df_pairs, lineage_A=lineage_A, lineage_B=lineage_B))
}

# function: build phased locus from the consensus and the block haplotypes (source: Claude)
f_build_phased_locus <- function(block_cds, df_blocks) {
    h0 <- rep("N", Biostrings::width(block_cds)[1])
    h1 <- h0

    # iterate over blocks
    for (i in 1:nrow(df_blocks)) {
        prefix   <- paste0(df_blocks$locus[i], "_PS", df_blocks$PS[i])

        # extract haplotype blocks
        block_h1 <- strsplit(toupper(as.character(block_cds[[paste0(prefix, "_h1")]])), "")[[1]]
        block_h2 <- strsplit(toupper(as.character(block_cds[[paste0(prefix, "_h2")]])), "")[[1]]
        covered  <- block_h1 != "N"

        # keep: block h1 is lineage A; flip: block h2 is lineage A
        if (df_blocks$decision[i] == "keep") {
            h0[covered] <- block_h1[covered]
            h1[covered] <- block_h2[covered]
        } else if (df_blocks$decision[i] == "flip") {
            h0[covered] <- block_h2[covered]
            h1[covered] <- block_h1[covered]
        }
    }

    return(Biostrings::BStringSet(c(h0=paste(h0, collapse=""), h1=paste(h1, collapse=""))))
}

# function: run AMAS
f_amas_dna <- function(input_regex, fn_output_aln, fn_output_part, exe_amas) {
    cmd_amas <- paste(exe_amas, "concat -f fasta -d dna",
                      "-i", input_regex,
                      "--part-format nexus",
                      "-t", fn_output_aln,
                      "-p", fn_output_part)
    system(cmd_amas)
}

# function: convert "N" and "?" to gaps
f_unknown2gap <- function(fn_fasta, fn_output) {
    # read the DNA alignment
    seq <- Biostrings::readBStringSet(fn_fasta, format="fasta")

    # convert "N" to gaps
    for (i in seq_along(seq)) {
        seq[i] <- gsub("N", "-", seq[i])
        seq[i] <- gsub("\\?", "-", seq[i])
    }

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)
}

# function: extract best models for each partition
f_extract_best_models <- function(fn_best_models, fn_partition) {
    # open nexus file
    ls_partition <- readLines(fn_partition)

    # open best model per locus
    df_model <- data.table::fread(fn_best_models)

    # extract the partition names from the nexus file
    ls_charset <- grep("charset", ls_partition, value=TRUE, ignore.case=TRUE)
    ls_charset <- lapply(ls_charset, function(line) {
        locus <- str_match(line, "charset\\s+(\\S+)\\s*=")[,2]
        coords <- str_match(line, "=\\s*(.+?);")[,2]
        data.frame(locus=locus, coords=str_trim(coords))
    })

    # combine the results into a single data frame
    df_charset <- do.call(rbind, ls_charset)

    # update locus name
    df_charset$locus <- sapply(df_charset$locus, function(x) {
        ls_parts <- unlist(strsplit(x, split="_"))
        ls_parts[length(ls_parts)]
    })

    # merge the two data.frames
    df_merge <- merge(df_charset, df_model, by="locus", all.x=TRUE)

    return(df_merge)
}

# function: run IQ-Tree 2 (fix topology)
f_iqtree2_fixtop <- function(fn_input, fn_tree, fn_partition, prefix, thread, is_redo, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-g", fn_tree,
                         "--prefix", prefix,
                         "-T", thread, "--quiet")

    if (!is.null(fn_partition) && fn_partition != "") {
        cmd_iqtree2 <- paste(cmd_iqtree2, "-p", fn_partition)
    }        

    if (is_redo) {
        cmd_iqtree2 <- paste(cmd_iqtree2, "-redo")
    }

    system(cmd_iqtree2)
}