# Example on meelup-mallee short reads
Meelup mallee (*E. x phylacis*) is an F1 hybrid between two eucalypt species, which are hypothesised to be *E. decipiens* and *E. virginea*. Here, we tested different parameters and assessed if we could consistenly recover the two parental lineages across loci. 

## Table of Content
- <a href="#01">Raw Data</a>
- <a href="#02">Running CAPTUS</a>
- <a href="#03">Running ShortReadsPhaser</a>
- <a href="#04">Results: Mapping</a>
- <a href="#05">Results: Species Tree Inference</a>
    - <a href="#05a">60x-coverage data</a>
    - <a href="#05b">30x-coverage data</a>

## <a id="01">Raw Data</a>
The raw short-read data are located at `/home/raw_data/MM/2018_short_reads/`. Each sample is labelled as `MMAA`, `MMAB`, `MMAC`, `MMBA`, etc., where:
- `MM` stands for meelup mallee
- The next letter is the branch ID (A-T)
- The final letter is the leaf ID (A-C)

Each leaf sample has ~6x coverage (based on an average genome size of 500Mb), which means that each branch has roughly 18x coverage across three leaf samples. In addition, each sample was sequenced twice, resulting in ~36x total coverage for individual branches. Among these branches, there are three stems where the sample size is doubled:
- Stem 17: `MMK` and `MMQ`
- Stem 18: `MML` and `MMR`
- Stem 23: `MMS` and `MMT`

Here, we combined the `MML` and `MMR` short reads using `cat` command in Unix, so we have two datasets with ~30x and ~60x coverage. Running `seqkit stats *.fastq.gz` generated the following:
- ~30x coverage: `SN877_0426_MML` + `SN8777_0429_MMR`
    ```
    file                      num_seqs         sum_len  min_len  avg_len  max_len
    phylacis_R1.fastq.gz    55,992,057   8,454,800,607      151      151      151
    phylacis_R2.fastq.gz    55,992,057   8,454,800,607      151      151      151
    ```
- ~60x coverage: `SN877_0426_MML` + `SN877_0427_MML` + `SN8777_0429_MMR` + `SN8777_0430_MMR`
    ```
    file                      num_seqs         sum_len  min_len  avg_len  max_len
    phylacis_R1.fastq.gz   125,404,800  18,936,124,800      151      151      151
    phylacis_R2.fastq.gz   125,404,800  18,936,124,800      151      151      151
    ```

## <a id="02">Running CAPTUS</a>
We ran <a href="https://doi.org/10.1093/sysbio/syag046">CAPTUS</a> on the paired-end reads to extract 1,187 BUSCO loci from <a href="https://doi.org/10.64898/2026.09.25.753384">Ivan et al. (2026)</a> with default parameters. We used the `-d` flag during `captus extract` to set the target loci as miscellaneous DNA markers because the BUSCO sequences are not codons.

> [!NOTE]
> Dataset from Ivan et al. (2026) does not have any *E. virginea* sample. The closest relative to *E. virginea* comes from *E. relicta*.

In brief, CAPTUS starts by trimming the adaptors and cleaning the reads using BBTools. Then, it *de novo* assembles the reads using MEGAHIT, and extracts the target loci from the assembly using BLAT. This step generates two output files that are relevant to this pipeline:
- `DNA_matches.fna`: individual hits of the target sequences
- `DNA_matches_flanked.fna`: individual hits of the target sequences + additional flanking sequences

## <a id="03">Running ShortReadsPhaser</a>
In `ShortReadsPhaser` pipeline, there are several free parameters to be set by users. Following the specification of our dataset, we used the following values as our `default` parameters (please refer to <a href="../config.yaml">`config.yaml`</a>):
```
# use individual gene trees inferred from Ivan et al. (2026) as backbone trees
dir_backbone_genetrees <- "/home/jeremias/captus/701_sample/tree_treefile/"

# use the CAPTUS output
tool <- "captus"

# set the minimum and maximum depths to call a variant
min_depth <- 15     # for 60x-coverage data, and 7.5 for 30x-coverage data
max_depth <- 150    # for 60x-coverage data, and 75 for 30x-coverage data

# in this case, the third rank represents series-level where *E. decipiens* and *E. virginea* differ
taxonomy_rank <- 3

# set the minimum accumulated likelihood weight ratio (alwr) to confidently assign a haplotype placement
min_alwr <- 0.8

# set the minimum fraction of placements that must agree to call a lineage
min_consensus_thresh <- 1.0
```

## <a id="05">Results: Mapping</a>

### Calculating expected pairwise genetic distance between *E. decipiens* and *E. relicta*
First, we calculated the pairwise genetic distance between `decipiens_CCA1310` (i.e., a representative sample of *E. decipiens*) and `relicta_CCA5717` (i.e., a closely-related species with *E. virginea*) using the BUSCO alignments from Ivan et al. (2026):
```
library(ape)

files  <- list.files("/home/jeremias/captus/busco_alignments", full.names=TRUE)
ls_ref <- c("decipiens_CCA1310", "relicta_CCA5717")

# iterate over BUSCO loci
res <- do.call(rbind, lapply(files, function(f) {
  # open the BUSCO alignment
  aln <- read.dna(f, format="fasta", as.character=TRUE)
  if (!all(ls_ref %in% rownames(aln))) {
    return(NULL)
  }

  # extract the two reference samples
  m <- aln[ls_ref, , drop=FALSE]
  m <- tolower(m)

  # count the number of sites with no missing data
  ok <- colSums(matrix(m %in% c("a", "c", "g", "t"), nrow=2)) == 2

  # return the number of sites per locus
  data.frame(locus = basename(f),
             diffs = sum(m[1, ok] != m[2, ok]),
             sites = sum(ok))
}))

# calculate the pairwise genetic distance per locus
res$p <- res$diffs / res$sites
```

Across 1,186 BUSCO where *both* were found, the average pairwise genetic distance between the two samples is 1.5% (SD=0.7%; Fig.1). This means that there is on average **one heterozygous site every ~70 bp**.

![Figure 1. Distribution of pairwise genetic distance between decipiens_CCA1310 and relicta_CCA5717 across 1,186 BUSCO loci.](./figures/fig2.png)
*<b>Figure 1.</b> Distribution of pairwise genetic distance between `decipiens_CCA1310` and `relicta_CCA5717` across 1,186 BUSCO loci. Vertical red line denotes the average across loci.*

### Mapping 60x-coverage data to the 1,187 CAPTUS best hits
After running CAPTUS on the 60x-coverage data, we mapped the cleaned short reads to the 1,187 CAPTUS best hits. We then assessed the mapping statistics as follows:

```
$ samtools idxstats phylacis.dedup.bam | awk '{T+=$2} END {print T}'
  6160664

$ samtools flagstat phylacis.sorted.bam
  249346330 + 0 in total (QC-passed reads + QC-failed reads)
  249094936 + 0 primary
  0 + 0 secondary
  251394 + 0 supplementary
  0 + 0 duplicates
  0 + 0 primary duplicates
  13305405 + 0 mapped (5.34% : N/A)
  13054011 + 0 primary mapped (5.24% : N/A)
  249094936 + 0 paired in sequencing
  124547468 + 0 read1
  124547468 + 0 read2
  8029416 + 0 properly paired (3.22% : N/A)
  8627550 + 0 with itself and mate mapped
  4426461 + 0 singletons (1.78% : N/A)
  571410 + 0 with mate mapped to a different chr
  446208 + 0 with mate mapped to a different chr (mapQ>=5)

$ samtools flagstat phylacis.dedup.bam
  6895181 + 0 in total (QC-passed reads + QC-failed reads)
  6895181 + 0 primary
  0 + 0 secondary
  0 + 0 supplementary
  0 + 0 duplicates
  0 + 0 primary duplicates
  6895181 + 0 mapped (100.00% : N/A)
  6895181 + 0 primary mapped (100.00% : N/A)
  6895181 + 0 paired in sequencing
  3457041 + 0 read1
  3438140 + 0 read2
  6382983 + 0 properly paired (92.57% : N/A)
  6733591 + 0 with itself and mate mapped
  161590 + 0 singletons (2.34% : N/A)
  329849 + 0 with mate mapped to a different chr
  329849 + 0 with mate mapped to a different chr (mapQ>=5)
```

GenomeScope2 estimated the genome size of meelup mallee to be roughly <a href="https://github.com/roblanf/MM_genome">530Mb</a> based on long sequencing reads, which means that the total length of the 1,187 BUSCO loci (~6.1Mb) accounts for ~1.2% of the genome. After deduplication, the retained reads represent ~2.8% of the total reads, which is ~2.3x higher than the expectation. In other words, the filtered BAM has roughly half off‑target reads (i.e., reads that do not come from the BUSCO loci). 

We then calculated the per-locus average mapping depth to assess if the extra mapped reads are localised or spread evenly across all loci. To do this, we first calculated the mapping depth for each CDS block using `samtools depth -a phylacis.dedup.bam`. Note that a single locus may contain multiple CDS blocks. The median mapping depth across all blocks is 62x, but the average is 105.89x, which reflects that the off-target reads are concentrated in a minority of positions. If we calculated the per-locus average mapping depth across all CDS blocks, there are 133 BUSCO loci with >=150 average depth (Fig. 2).

![Figure 2. Per-locus average mapping depth of the meelup mallee short-read data.](./figures/fig1.png)
*<b>Figure 2.</b> Per-locus average mapping depth of the meelup mallee short-read data. Vertical red dashed lines denote the minimum and maximum depths during variant calling.*

## <a id="05">Results: Species Tree Inference</a>

### <a id="05a">60x-coverage data</a>
We first retrieved all CDS blocks across the 1,187 BUSCO loci, and filtered out four loci with all Ns after remapping (including `167757at71240` with no haplotype recovered):

| Step         | #BUSCO | #Blocks |
| ------------ | -----: | ------: |
| Initial      | 1,187  | 4,713   |
| Filtered out | 4      | 7       |  
| Passed       | 1,183  | 4,676   |

After running EPA-NG, we found that the most common parental lineages recovered from the gene trees are `Curviptera` and `Falcatae`. Across the 4,676 blocks, there were 3,205 blocks unresolved at series-level (`taxonomic_rank==3`), while the rest (1,471 blocks) can be categorised as:

| Criterion                                       | #Blocks |
| ----------------------------------------------- | ------: |
| Either or both haplotypes have alwr >= min_alwr | 721     |
| Both haplotypes have alwr < min_alwr            | 743     | 
| Both haplotypes have the same lineages          | 7       |

Among 721 blocks, 117 blocks (across 114 BUSCO loci) have confident placement for both haplotypes. Thus, we inferred concatenated species trees using two sets of loci:
- **114 BUSCO** with confident placement for both haplotypes
- **671 BUSCO** with confident placement for either or both haplotypes

Both sets resulted in identical placement of the two haplotypes: one as sister to `decipiens_CCA1310` and the other as sister to `relicta_CCA5717`, while the unphased haplotype becomes the sister to `relicta_CCA5717` with a slightly longer branch length.
- Phylogenetic placements on `114-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.0036)
    - Haplotype two: sister to `relicta_CCA5717` (0.0064)
    - Unphased (CAPTUS best hits): sister to `relicta_CCA5717` (0.0093)
- Phylogenetic placements on `671-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.0044)
    - Haplotype two: sister to `relicta_CCA5717` (0.0084)
    - Unphased (CAPTUS best hits): sister to `relicta_CCA5717` (0.0090)

### <a id="05b">30x-coverage data</a>
We did the same analyses on the dataset with ~30x coverage, and recovered similar results:
- Number of BUSCO with confident placement for both haplotypes: **107**
- Number of BUSCO with confident placement for either or both haplotypes: **675**
- Phylogenetic placements on `107-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.0034)
    - Haplotype two: sister to `relicta_CCA5717` (0.0065)
    - Unphased (CAPTUS best hits): sister to `relicta_CCA5717` (0.0084)
- Phylogenetic placements on `675-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.0046)
    - Haplotype two: sister to `relicta_CCA5717` (0.0086)
    - Unphased (CAPTUS best hits): sister to `relicta_CCA5717` (0.0092)

---
*Last update: 19 October 2026 by Jeremias Ivan*