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
                         as.character(Biostrings::reverseComplement(Biostrings::BStringSet(ref))))

        # align the CDS to the flanked reference
        aln <- Biostrings::pairwiseAlignment(Biostrings::BStringSet(cds), Biostrings::BStringSet(query))

        # extract the aligned positions and calculate identity
        p <- unlist(strsplit(as.character(Biostrings::pattern(aln)), ""))     # aligned CDS
        q <- unlist(strsplit(as.character(Biostrings::subject(aln)), ""))     # aligned flanked

        # extract the starting position of the aligned sites
        ref_pos   <- cumsum(q != "-") + Biostrings::start(Biostrings::subject(aln)) - 1
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
        seq <- as.character(Biostrings::reverseComplement(Biostrings::BStringSet(seq)))
    }

    # extract the coding sequence
    cds <- Biostrings::BStringSet(setNames(paste(substring(seq, exons$start, exons$end), collapse=""), names(hap)[1]))
    return(cds)
}