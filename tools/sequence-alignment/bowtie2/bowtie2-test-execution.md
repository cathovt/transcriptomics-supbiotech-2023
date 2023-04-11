The most known algorithms for mapping are **bowtie2**,  **bwa** and **minimap**. There are also some others aligners more specific to RNA-seq such as **STAR** or **salmon**. For bacteria where the genomes are small with no repeated regions, 
those are not necessary.

In this work, **bowtie2** is used as **bwa** and **minimap** are not well suited for very short reads (< 100bp).

### Building indexes

The aligner algorithm needs to index fasta files. An **index** is similar to a dictionary in which the sequences present in the genome are sorted. 
Indexing allows to find the matching sequences very fast. 

To build the index using bowtie2:
> bowtie2-build 

To simplify the analysis we will just look at the genome of the bacteria at T0.

> cd transcriptomics-supbiotech-2023
> mkdir -p reference/index
> bowtie2-build reference/YB886.fa reference/index/YB886

Once the index is created, it is possible to map the reads with --very-sensitive method [Documentation](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
> bowtie2 --help

### Seed

As the simultaneous mapping of all sequence is usually slow for the aligner, it uses a subset of the sequence called **seed** to find the first potential position. Afterwards, the aligner extends the seed to find the best match. 

The seed is set up by the **-L parameter** in bowtie2. The bigger the seed is, the faster the algorithm is. Indeed, the aligner needs to extend on fewer positions, but it may miss the perfect position if there is no appropriate sequence in the seed.

### Alignment modes

**Paired-end alignment** typically means keeping track and reporting the alignment of the booth pairs in a read pair. Each read is aligned separately, and the information on both pairs is combined and reported in the same alignment line.

In **end-to-end mode**, an entire read must be aligned to the reference genome for an alignment to be returned. End-to-end alignments are considered a stricter method of alignment compared to local alignments.

Opposed to **end-to-end mode**, **local mode** allows the aligner to avoid mapping of the reads' edges.
It is usually advised to use it for RNA-seq as there is some splicing introducing issues with local alignment. 
However, as splicing does not occur in bacteria, **it is not necessary to use it in this work**

> mkdir -p results/alignment
> bowtie2 --very-sensitive --maxins 1000 -x reference/index/YB886 -1 results/cutadapt/Bacillus_infection_T0_R1_trimmed.fq.gz -2 results/cutadapt/Bacillus_infection_T0_R2_trimmed.fq.gz -S results/alignment/T0.sam



