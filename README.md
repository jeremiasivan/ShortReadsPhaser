# ShortReadsPhaser

**Short Reads Phaser** is an R pipeline to phase individual loci using short-read data. It requires CAPTUS output from the new sample, as well as backbone locus alignments and trees as reference. It is mainly developed and tested using Linux, so there might be incompatibilities using MacOS and Windows.

## Table of Content
- <a href="#prereqs">Prerequisites</a>
- <a href="#genpipe">General Pipeline</a>
- <a href="#refs">References</a>

## <a id="prereqs">Prerequisites</a>
Short Reads Phaser requires a number of software and R packages to run. We recommend you to use environment management system (e.g. `conda`) to install the prerequisites, but you can also use `install.packages()` built-in function in R or RStudio.

### Software
|    Name        |                             Website                              |                             Anaconda                             |
| -------------- |:----------------------------------------------------------------:|:----------------------------------------------------------------:|
| AMAS           | <a href="https://github.com/marekborowiec/AMAS">Link</a>         | <a href="https://anaconda.org/bioconda/amas">Link</a>            |
| BCFtools       | <a href="https://github.com/samtools/bcftools">Link</a>          | <a href="https://anaconda.org/bioconda/bcftools">Link</a>        |
| BWA            | <a href="https://github.com/lh3/BWA">Link</a>                    | <a href="https://anaconda.org/bioconda/bwa">Link</a>             |
| EPA-NG         | <a href="https://github.com/pierrebarbera/epa-ng">Link</a>       | <a href="https://anaconda.org/bioconda/epa-ng">Link</a>          |
| Gappa          | <a href="https://github.com/lczech/gappa">Link</a>               | <a href="https://anaconda.org/bioconda/gappa">Link</a>           |
| IQ-TREE        | <a href="http://www.iqtree.org">Link</a>                         | <a href="https://anaconda.org/bioconda/iqtree">Link</a>          |
| MAFFT          | <a href="https://github.com/GSLBiotech/mafft">Link</a>           | <a href="https://anaconda.org/bioconda/mafft">Link</a>           |
| SAMtools       | <a href="https://github.com/samtools/samtools">Link</a>          | <a href="https://anaconda.org/bioconda/samtools">Link</a>        |
| WhatsHap       | <a href="https://github.com/whatshap/whatshap">Link</a>          | <a href="https://anaconda.org/bioconda/whatshap">Link</a>        |

### R packages
|    Name    |                               CRAN                               |                                   Anaconda                               |
| ---------- |:----------------------------------------------------------------:|:------------------------------------------------------------------------:|
| Biostrings | <a href="https://bioconductor.org/packages/Biostrings">Link</a>  | <a href="https://anaconda.org/bioconda/bioconductor-biostrings">Link</a> |
| data.table | <a href="https://cran.r-project.org/package=data.table">Link</a> | <a href="https://anaconda.org/conda-forge/r-data.table">Link</a>         |
| doSNOW     | <a href="https://cran.r-project.org/package=doSNOW">Link</a>     | <a href="https://anaconda.org/conda-forge/r-dosnow">Link</a>             |
| log4r      | <a href="https://cran.r-project.org/package=log4r">Link</a>      | <a href="https://anaconda.org/conda-forge/r-log4r">Link</a>              |
| optparse   | <a href="https://cran.r-project.org/package=optparse">Link</a>   | <a href="https://anaconda.org/conda-forge/r-optparse">Link</a>           |
| rmarkdown  | <a href="https://cran.r-project.org/package=rmarkdown">Link</a>  | <a href="https://anaconda.org/conda-forge/r-rmarkdown">Link</a>          |
| tidyverse  | <a href="https://cran.r-project.org/package=tidyverse">Link</a>  | <a href="https://anaconda.org/conda-forge/r-tidyverse">Link</a>          |
| yaml       | <a href="https://cran.r-project.org/package=yaml">Link</a>       | <a href="https://anaconda.org/conda-forge/r-yaml">Link</a>    

## <a id="genpipe">General Pipeline</a>
1. **Clone the Git repository** <br>
    ```
    git clone git@github.com:jeremiasivan/ShortReadsPhaser.git
    ```

2. **Install the prerequisites** <br>
    - Create a new conda environment
        ```
        conda create -n taxontest
        conda activate taxontest
        ```
    - Installing prerequisites
        ```
        conda install -c conda-forge r-data.table r-doSNOW r-log4r r-optparse r-rmarkdown r-tidyverse r-yaml
        ```

3. **Update the parameters in `config.yaml`** <br>

4. **Run Short Reads Phaser** <br>
    ```
    Rscript run_pipeline.R --config config.yaml
    Rscript run_pipeline.R --config config.yaml --redo
    ```

    In UNIX-based operating systems (e.g., Linux and MacOS), it is advisable to use `nohup` or `tmux` to run the whole pipeline. For Windows, you can use `psmux`. 

---
## <a id="refs">References</a>
1. Borowiec, M.L. (<a href="https://doi.org/10.7717/peerj.1660">2016</a>). **AMAS: A fast tool for alignment manipulation and computing of summary statistics**. *PeerJ*, *4*, e1660.

2. Danecek, P., et al. (<a href="https://doi.org/10.1093/gigascience/giab008">2021</a>). **Twelve years of SAMtools and BCFtools**. *GigaScience*, *10*(2), giab008.

3. Li, H. & Durbin, R. (<a href="https://doi.org/10.1093/bioinformatics/btp324">2009</a>). **Fast and accurate short read alignment with Burrows–Wheeler transform**. *Bioinformatics*, *25*(14), 1754–1760.

4. Barbera, P., et al. (<a href="https://doi.org/10.1093/sysbio/syy054">2019</a>). **EPA-ng: Massively Parallel Evolutionary Placement of Genetic Sequences**. *Systematic Biology*, *68*(2), 365–369.

5. Czech, L., et al. (<a href="https://doi.org/10.1093/bioinformatics/btaa070">2020</a>). **Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data**. *Bioinformatics*, *36*(10), 3263–3265.

6. Minh, B.Q., et al. (<a href="https://doi.org/10.1093/molbev/msaa015">2020</a>). **IQ-TREE 2: New models and efficient methods for phylogenetic inference in the genomic era**. *Molecular Biology and Evolution*, *37*(5), 1530–1534.

7. Katoh, K. & Standley, D.M. (<a href="https://doi.org/10.1093/molbev/mst010">2013</a>). **MAFFT multiple sequence alignment software version 7: Improvements in performance and usability**. *Molecular Biology and Evolution*, *30*(4), 772–780.

8. Martin, M., et al. (<a href="https://doi.org/10.1101/085050">2016</a>). **WhatsHap: Fast and accurate read-based phasing**. *bioRxiv*.

9. Pagès, H., et al. (<a href="https://doi.org/10.18129/B9.bioc.Biostrings">2026</a>). **Biostrings: Efficient manipulation of biological strings**. *R package*.

10. Barrett, T., et al. (<a href="https://doi.org/10.32614/CRAN.package.data.table">2026</a>). **data.table: Extension of 'data.frame'**. *R package*.

11. Daniel, F. (<a href="https://cran.r-project.org/package=doSNOW">2022</a>). **doSNOW: Foreach Parallel Adaptor for the 'snow' Package**. *R package*.

12. White, J.M. & Jacobs, A. (<a href="https://doi.org/10.32614/CRAN.package.log4r">2024</a>). **log4r: A Fast and Lightweight Logging System for R, Based on 'log4j'**. *R package*.

13. Davis, T.L. (<a href="https://doi.org/10.32614/CRAN.package.optparse">2026</a>). **optparse: Command Line Option Parser**. *R package*.

14. Allaire, J.J., et al. (<a href="https://doi.org/10.32614/CRAN.package.rmarkdown">2026</a>). **rmarkdown: Dynamic Documents for R**. *R package*.

15. Wickham, H., et al. (<a href="https://doi.org/10.21105/joss.01686">2019</a>). **Welcome to the tidyverse**. *Journal of Open Source Software*, *4*(43), 1686.

16. Stephens, J., et al. (<a href="https://doi.org/10.32614/CRAN.package.yaml">2025</a>). **yaml: Methods to Convert R Data to YAML and Back**. *R package*.

17. Anthropic. (<a href="https://claude.ai/">2026</a>). Claude 5.1 Fable was used to generate `config.yaml` and `run_pipeline.R`. 

---
*Last update: 16 October 2026 by Jeremias Ivan*