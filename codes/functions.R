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
f_locate_exons <- function(fn_cds, fn_flanked) {
    # open the sequences
    cds <- toupper(as.character(Biostrings::readBStringSet(fn_cds)[[1]]))
    ref <- toupper(as.character(Biostrings::readBStringSet(fn_flanked)[[1]]))

    # check if flanked sequence is + or - strand
    best <- NULL
    for (strand in c("+", "-")) {
        # check the strand and reverse complement if necessary
        query <- ifelse (strand == "+",
                         ref,
                         as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(ref))))

        # align the CDS to the flanked reference
        aln <- pwalign::pairwiseAlignment(Biostrings::BStringSet(cds), Biostrings::BStringSet(query), type="global-local")

        # extract the aligned positions and calculate identity
        p <- unlist(strsplit(as.character(pwalign::pattern(aln)), ""))     # aligned CDS
        q <- unlist(strsplit(as.character(pwalign::subject(aln)), ""))     # aligned flanked

        # extract the starting position of the aligned sites
        ref_pos   <- cumsum(q != "-") + pwalign::start(pwalign::subject(aln)) - 1
        mapped    <- p != "-" & q != "-"
        positions <- ref_pos[mapped]
        identity  <- mean(p[mapped] == q[mapped])
        if (is.null(best) || length(positions) > length(best$positions)) {
            best <- list(strand=strand, positions=positions, identity=identity, unmapped=sum(p != "-" & q == "-"))
        }
    }

    # accept only if high identity
    if (length(best$positions) < 0.9*nchar(cds) || best$identity < 0.95) {
        return(NULL)
    }

    # extract consecutive positions
    breaks <- which(diff(best$positions) != 1)
    df_output <- data.table::data.table(strand=best$strand, identity=round(best$identity, 4), unmapped=best$unmapped,
                                        start=best$positions[c(1, breaks + 1)], end=best$positions[c(breaks, length(best$positions))])

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
                       "--per-query-results --best-hit --allow-file-overwriting",
                       "--out-dir", outdir,
                       ">>", log_file)
    system(cmd_gappa)
}

# function: read gappa assignments
f_read_gappa_assignments <- function(dir_epa_ng, rank) {
    # list output files from gappa examine assign
    fn_per_query <- list.files(dir_epa_ng, pattern="^per_query\\.tsv$", recursive=TRUE, full.names=TRUE)

    # combine all assignments into a data.frame
    ls_assign <- lapply(fn_per_query, function(file) {
        df <- data.table::fread(file)
        df$locus <- basename(dirname(dirname(file)))
        df
    })
    df_assign <- data.table::rbindlist(ls_assign)

    # extract information from the sequence name
    df_assign$PS      <- gsub(".*_PS([^_]+)_h[12]$", "\\1", df_assign$name)
    df_assign$hap     <- gsub(".*_(h[12])$", "\\1", df_assign$name)
    df_assign$lineage <- sapply(strsplit(df_assign$taxopath, split=";"), function(x) { x[min(rank, length(x))] })
    df_assign$aLWR    <- as.numeric(df_assign$aLWR)

    # subset the columns of interest
    df_assign <- df_assign %>% select(c("locus", "PS", "hap", "lineage", "aLWR"))

    return(df_assign)
}

# function: find the two parental lineages for every block (source: Claude)
f_assign_blocks <- function(df_assign, min_alwr) {
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
    lineage_A <- most_common_pair[1]
    lineage_B <- most_common_pair[2]

    # assign decisions per block
    df_blocks$decision <- "ambiguous"
    for (i in 1:nrow(df_blocks)) {
        # extract the lineages of the two haplotypes
        l1 <- df_blocks$h1_lineage[i]
        l2 <- df_blocks$h2_lineage[i]

        # assign decisions based on the lineages and confidence
        if (is_both[i] && l1 == lineage_A && l2 == lineage_B) {
            df_blocks$decision[i] <- "keep"
        } else if (is_both[i] && l1 == lineage_B && l2 == lineage_A) {
            df_blocks$decision[i] <- "flip"
        } else if (is_both[i] && l1 == l2) {
            df_blocks$decision[i] <- "same"
        } else if (is_conf1[i] && l1 %in% c(lineage_A, lineage_B) && !(l2 %in% c(lineage_A, lineage_B))) {
            df_blocks$decision[i] <- ifelse(l1 == lineage_A, "keep", "flip")      
        } else if (is_conf2[i] && l2 %in% c(lineage_A, lineage_B) && !(l1 %in% c(lineage_A, lineage_B))) {
            df_blocks$decision[i] <- ifelse(l2 == lineage_A, "flip", "keep")      
        }
    }

    return(list(blocks=df_blocks, pairs=df_pairs, lineage_A=lineage_A, lineage_B=lineage_B))
}