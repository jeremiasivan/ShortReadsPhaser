# Example on meelup-mallee short reads
Meelup mallee (*E. x phylacis*) is an F1 hybrid between two eucalypt species, that are hypothesised to be *E. decipiens* and *E. virginea*. Here, we tested different parameters and assessed if we could consistenly recover the two parental lineages across loci. 

## Table of Content
- <a href="#01">Raw Data</a>
- <a href="#02">Running CAPTUS</a>
- <a href="#03">Running ShortReadsPhaser</a>
- <a href="#04">Results</a>
- <a href="#05">Using Other Parameters</a>

## <a id="01">Raw Data</a>
The raw short-read data are located at `/home/raw_data/MM/2018_short_reads/`. Each sample is labelled as `MMAA`, `MMAB`, `MMAC`, `MMBA`, etc., where:
- `MM` stands for meelup mallee
- The next letter is the branch ID (A-T)
- The final letter is the leaf ID (A-C)

Each leaf sample has ~12x coverage (based on an average genome size of 500Mb), which means that each branch has roughly 36x coverage across three leaf samples. Among these branches, there are three stems where the sample size is doubled:
- `MMK` and `MMQ`
- `MML` and `MMR`
- `MMS` and `MMT`

Here, we combined the `MML` and `MMR` short reads using `cat` command in Unix, resulting in paired-end reads of ~60x coverage.

## <a id="02">Running CAPTUS</a>
We ran <a href="https://doi.org/10.1093/sysbio/syag046">CAPTUS</a> on the paired-end reads to extract 1,187 BUSCO loci from <a href="https://doi.org/10.64898/2026.09.25.753384">Ivan et al. (2026)</a> with default parameters. We used the `-d` flag during `captus extract` to set the target loci as miscellaneous DNA markers because the BUSCO sequences do not follow the codon positions.

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

# set the minimum and maximum depths to call a variant based on ~60x average coverage
min_depth <- 15
max_depth <- 160

# in this case, the third rank represents series-level where *E. decipiens* and *E. virginea* differ
taxonomy_rank <- 3

# set the minimum accumulated likelihood weight ratio (alwr) to confidently assign a haplotype placement
min_alwr <- 0.8
```

## <a id="04">Results</a>
We first retrieved all CDS blocks across the 1,187 BUSCO loci, and excluded those with all Ns after remapping.

| Step     | #BUSCO | #Blocks |
| -------- | -----: | ------: |
| Initial  | 1,187  | 4,684   |
| All Ns   | 3      | 8       |  
| Passed   | 1,184  | 4,676   |

After running EPA-NG, we found that the most common parental lineages recovered from the gene trees are `Curviptera` and `Falcatae`. Across the 4,676 blocks:

| Criterion                                       | #Blocks |
| ----------------------------------------------- | ------: |
| Either or both haplotypes have alwr >= min_alwr | 735     |
| Both haplotypes have alwr < min_alwr            | 827     | 
| Both haplotypes have the same lineages          | 5       |

Among 735 blocks, 106 blocks (across 102 BUSCO loci) have confident placement for both haplotypes. Thus, we inferred concatenated species trees using two sets of loci:
- **102 BUSCO** with confident placement for both haplotypes
- **663 BUSCO** with confident placement for either or both haplotypes

Both sets resulted in identical placement of the two haplotypes: one as sister to `decipiens_CCA1310` and the other as sister to `relicta_CCA5717`. The `102-BUSCO` tree has slightly shorter branch lengths for both haplotypes (0.0043 and 0.0065, respectively) compared to the `663-BUSCO` tree (0.0052 and 0.0085, respectively).

## <a id="05">Using Other Parameters</a>

> [!NOTE]
> ✅: consistent output using `default` parameters
> ❌: inconsistent with output using `default` parameters

### `min_alwr = 0.7`
- Number of BUSCO with confident placement for both haplotypes: **122**
- Number of BUSCO with confident placement for either or both haplotypes: **700**
- Phylogenetic placements on `122-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.004) ✅
    - Haplotype two: sister to `relicta_CCA5717` (0.0065) ✅
- Phylogenetic placements on `700-BUSCO` tree (branch lengths):
    - Haplotype one: sister to `decipiens_CCA1310` (0.0052) ✅
    - Haplotype two: sister to `relicta_CCA5717` (0.0085) ✅
