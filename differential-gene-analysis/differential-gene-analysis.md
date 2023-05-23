# Differential gene analysis

### Objective
Analyse RNAseq data and gene differential expression. The purpose of this analysis is to evaluate the expression of genes and to compare it across several samples.

### Data and requirements
Data from the RNAseq libraries of _Bacillus subtilis_ infected by its phage SPP1. 

**Reference genomes:** YB886 strain from _Bacillus subtilis_ and its phage SPP1. 

In this example, the **featureCount** table has been computed based on all the reads. 

There are 3 replicates for each of the 4 time points: 0, 5, 10, 25 minutes.

This section provides an R script for the analysis.

### Installation of R packages

R is installed from: http://cran.r-project.org/

In PyCharm, `R Language for IntelliJ` plugin enables processing of R files and functioning of packages.

In R console:

    install.packages("dplyr")
    install.packages("ggplot2")
    install.packages("ggbeeswarm")
    install.packages("ashr")
    install.packages("hwriter")
    install.packages("ini")
    install.packages("ggrepel")
    install.packages("png")
    install.packages("circlize")
    install.packages('BiocManager')
    BiocManager::install("DESeq2")
    BiocManager::install("tximport")
    BiocManager::install("ReportingTools")
    BiocManager::install("ComplexHeatmap") 
    BiocManager::install("rhdf5")

In addition, [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889) software is used.

To install, download the package from [SourceForge] (https://sourceforge.net/projects/subread/files/).

Unpack the archive with `tar zxvf subread-2.x.x.tar.gz`, where `.x.x` is the version of subRead. 
Copy the binary executables from `bin` folder to `/usr/local/bin`. 

In terminal, run:

    featureCounts --help

## Computing gene counts

The first step for a differential gene expression analysis is to **count the numbers of reads/fragments mapping on a gene**. 

To do that, the positions of genes are required. 

Usually, genome annotation files with the coding sequence position are provided. 
Sometimes it is necessary to generate them by yourself, for example, using [prokka](https://github.com/tseemann/prokka).
In this case,As _Bacillus subtilis_ have been extensively annotated, it yields a satisfying annotation.

### .gff 

The annotations are written in the annotation file. 
Here, we use a `.gff` file format, however several other formats exist. The main ones are: `.gff`, `.gtf` and `.genbank`. 

`.gff` and `.gtf` are quite similar with 1 entry per annotation, while `.genbank` has several lines for each annotation starting with a field.

In the command line:

    cd # Directory containing the .gff file
    head -n 20 reference/YB886.gff 

The `.gff` file is composed of **a header with lines starting with ##**. 
In the header one will usually find the **version of the .gff format**, the **sequence names** and **length of the original fasta**. 

Afterwards, there are annotations with 1 line per each. Each line has 9 columns as described in the table below:
![Снимок экрана 2023-05-22 в 10.01.53.png](..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_QI32zd%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-05-22%20%D0%B2%2010.01.53.png)

The last one is a column where one can add any additional field. 
Some pipelines will need a specific field, such as a `gene_id field` for `featureCounts`.

### Running featureCounts

    cd Desktop/transcriptomics_supbiotech_2023/results
    mkdir -p featureCounts
    featureCounts \
    -p \
    -O \
    -s 1 \
    -t CDS \
    --ignoreDup \
    -T 1 \
    -G reference/YB886.fa \
    -a reference/YB886.gff \
    -o results/featureCounts/tmp_counts.txt \
    alignment/compressed-results/T0_sorted.bam alignment/compressed-results/T10_sorted.bam

    # To format the header:

    FILE_PATH=results/alignment/
    sed "s|$FILE_PATH||g" results/featureCounts/tmp_counts.txt
    sed 's/_sorted.bam//g' results/featureCounts/tmp_counts.txt > results/featureCounts/featureCounts.txt
    rm results/featureCounts/tmp_counts.txt

**Question 1:** what does the `-O` option mean?

If there is a read overlap between 2 (or more) different features (CDS), it will be counted twice (or more). 
As in bacteria genes could be organized in operons, the same RNA can overlap multiple genes. If it happens, we are interested to count it for both genes.

**Question 2:** what does the `-s 1` option mean ?

This options means that the library is stranded. 
During the RNA sequencing, the information of the number of RNA strands is kept. 
By giving this option, we ask only to count those reads which belong to the **1st** strand of the genome.

To look at the output file:

        head -n 6 results/featureCounts/featureCounts.txt

The values are quite low here as the sequence comes from the subsampled libraries.

## Differential gene expression analysis

### Running DESeq2

**Inputs:** **gene count table** and **a metadata file describing all the conditions**.
    
    transcriptomics_supbiotech_2023/metadata.tsv
    featureCounts_all.txt

## Results

### PCA plot

> PCA stands for **Principle Component Analysis**. It’s a dimensionality-reduction method, transforming a large dataset of many dimensions into a smaller dataset of fewer dimensions. 
> The idea is to look at the values of each dimension: 
> 
> - gene expression values 
> - variance between the samples 
> 
> and compute new dimensions which are combinations of the original ones. 

Here is the PCA plot of our sample:

![PCAplot.png](deseq2-results%2FPCAplot.png)

- **There is more variance between conditions than between replicates:** second component (PC2) demonstrates less variance. 
This is reassuring, proving that it is possible to see the actual difference between samples, and that the phage infection impacts the bacteria.

-  Second PCA component (PC2) demonstrates the noise between replicates, suggesting that they could not be perfectly reliable.

### MA plot

> **MA** stands for **Mean Average** plot. 
> X-axis refers to the **mean average count** while y-axis shows the **log fold change** between 2 experiments. 
> 
> The MA plot is a quality control: the dots are colored if they are considered as significant with **p-value below a given threshold**. 
> 
> The **p-value** depends on the variance between 2 conditions and the one within each condition. Thus, a high p-value (grey dot) can be due to either no variance between 2 conditions, or a high variance within one condition. 
> 
> For the perfect replicates, there shouldn’t be any overlap between blue and grey dots.

![time_0_vs_25_MAplot.png](deseq2-results%2Ftime_0_vs_25_MAplot.png)

- As there is no huge overlap between significant and non-significant values, the replicates seem fine for both conditions

- The grey area is larger at low mean of normalized counts due **to higher noise for low signals**. To be considered as significant, the value requires a **bigger fold change to overcome the noise variance**.

### The Volcano plot

> The volcano plot is a plot of the fold change depending on the p-value. The name of this plot originates from its shape. 
> 
> The plot aims at **separating differentially expressed genes from the others**. 
> 
> 2 thresholds are used: a **p-value threshold** specifying the log fold change that is relevant and significantly different from 0, and the **fold change**.

> **The upper left panel** represents the down-regulated genes with negative fold change.
> 
> **The upper right panel** represents up-regulated genes.

![time_0_vs_25_VolcanoPlot.png](deseq2-results%2Ftime_0_vs_25_VolcanoPlot.png)

**The meaning of the colours:**

- **Grey:** Non-significant fold change. Normal genes with no differential expression.
- **Yellow:** Huge fold change, but non-significant. Genes with high noise between replicates and conditions.
- **Red:** Significant, but small fold change. Significantly differentially expressed genes, but with only small variation. They are usually not kept due to possibly little biological relevance due to small variations.
- **Blue:** Huge and significant fold change. Significantly differentially expressed genes with huge variation, which are relevant.

![replicate_1_vs_2_VolcanoPlot.png](deseq2-results%2Freplicate_1_vs_2_VolcanoPlot.png)

In the volcano plot between 2 replicates, no genes are differentially expressed.

### Heatmap

> The heatmap is a plot representing **the counts of significantly differentially expressed genes across the conditions**. 
> It allows to see if there are different populations of genes among the significant ones.

![time_0_vs_5_heatmap.png](deseq2-results%2Ftime_0_vs_5_heatmap.png)

3 populations:

- Genes down-regulated across time
- Genes up-regulated across time
- Genes up-regulated at 5 minutes and down-regulated later.

# What next?

The next aim is to check the relevance and the functions of identified differentially expressed genes in the biological experiments.

To do so, it's possible to perform:

- Gene set enrichment analysis to look for differentially expressed pathways
- Bibliography research to look if the differentially expressed genes could have been expected
- Functions research for hypothetical proteins, using folding prediction (Alphafold), search for homology 
- Mutation of genes in the biological systems (bacteria) to see the impact
- Comparison with other genomic tracks if some are available

    