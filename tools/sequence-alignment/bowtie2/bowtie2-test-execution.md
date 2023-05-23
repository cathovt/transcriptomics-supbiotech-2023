The most known algorithms for mapping are **bowtie2**,  **bwa** and **minimap**. There are also some others aligners more specific to RNA-seq such as **STAR** or **salmon**. For bacteria where the genomes are small with no repeated regions, 
those are not necessary.

In this work, **bowtie2** is used as **bwa** and **minimap** are not well suited for very short reads (< 100bp).

### Building indexes

The aligner algorithm needs to index fasta files. An **index** is similar to a dictionary in which the sequences present in the genome are sorted. 
Indexing allows to find the matching sequences very fast. 

To build the index using bowtie2:
`bowtie2-build`

To simplify the analysis we will just look at the genome of the bacteria at T0.

>`cd transcriptomics-supbiotech-2023`
> 
>`mkdir -p reference/index`
> 
>`bowtie2-build reference/YB886.fa reference/index/YB886`

Once the index is created, it is possible to map the reads with `--very-sensitive` method. 

[Documentation:](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
`bowtie2 --help`

### Seed

As the simultaneous mapping of all sequence is usually slow for the aligner, it uses a subset of the sequence called **seed** to find the first potential position. Afterwards, the aligner extends the seed to find the best match. 

The seed is set up by the **-L parameter** in bowtie2. The bigger the seed is, the faster the algorithm is. Indeed, the aligner needs to extend on fewer positions, but it may miss the perfect position if there is no appropriate sequence in the seed.

### Alignment modes

> **Paired-end alignment** typically means keeping track and reporting the alignment of the booth pairs in a read pair. Each read is aligned separately, and the information on both pairs is combined and reported in the same alignment line.

In **end-to-end mode**, an entire read must be aligned to the reference genome for an alignment to be returned. End-to-end alignments are considered a stricter method of alignment compared to local alignments.

Opposed to **end-to-end mode**, **local mode** allows the aligner to avoid mapping of the reads' edges.
It is usually advised to use it for RNA-seq as there is some splicing introducing issues with local alignment. 
However, as splicing does not occur in bacteria, **it is not necessary to use it in this work**

### Alignment of test sequences

**First time-point (T0)**

> `mkdir -p alignment/results`
> 
> `bowtie2 --very-sensitive --maxins 1000 -x reference/index/YB886 -1 cutadapt/results/cutadapt/Bacillus_infection_T0_R1_trimmed.fq.gz -2 cutadapt/results/cutadapt/Bacillus_infection_T0_R2_trimmed.fq.gz -S alignment/results/T0.sam`

![Снимок экрана 2023-04-11 в 09.06.05.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_BSXMsZ%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-11%20%D0%B2%2009.06.05.png)

**Second time-point (T10)**

> `bowtie2 --very-sensitive --maxins 1000 -x reference/index/YB886 -1 cutadapt/results/cutadapt/Bacillus_infection_T10_R1_trimmed.fq.gz -2 cutadapt/results/cutadapt/Bacillus_infection_T10_R2_trimmed.fq.gz -S alignment/results/T10.sam`

![Снимок экрана 2023-04-11 в 09.14.14.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_foraGZ%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-11%20%D0%B2%2009.14.14.png)

**_Why the alignment rate is not 100%?_**
- There is a contamination with other DNA (another organism, vectors, adapters…). As tested sequences were checked in fastqc, the contamination is less likely.

- There are  sequencing errors which have introduced wrong sequences in the files, which now can’t be mapped. 

- There are mistakes in the reference genome, or it’s incomplete.

### The SAM format output

To show lines that don’t have @ in the beginning, to skip the header:

`grep -v '^@' alignment/results/T0.sam | head -n 4`

The file is a sam file, which is a tab separated file of 11 columns:
![Снимок экрана 2023-04-11 в 09.22.59.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_txLf4Y%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-11%20%D0%B2%2009.22.59.png)

The column that interest us most is the position of the reads on the reference (`RNAME`), the `FLAG`, information about the mapping and the mapping quality `MAPQ`, and the extent of sequence's uniqueness in the reference.
To decode a flag, we can use: [Decoding SAM Flags](https://broadinstitute.github.io/picard/explain-flags.html)

> **77**: read paired (0x1), read unmapped (0x4), mate unmapped (0x8), first in pair (0x40)
> 
> **141**: read paired (0x1), read unmapped (0x4), mate unmapped (0x8), second in pair (0x80)
> 
> **83**:  read paired (0x1), read mapped in proper pair (0x2), read reverse strand (0x10), first in pair (0x40)
> 
> **163**: read paired (0x1), read mapped in proper pair (0x2), mate reverse strand (0x20), second in pair (0x80)
> 
> **99**: read paired (0x1), read mapped in proper pair (0x2), mate reverse strand (0x20), first in pair (0x40)
> 
> **147**:    read paired (0x1), read mapped in proper pair (0x2), read reverse strand (0x10), second in pair (0x80)

`sam` files are quite heavy. To reduce their size, we could filter the mapping reads, sort them and index them using `samtools`. The  `bam` file is the compressed version of `sam`.

> - Sort and compress the `sam` file:
> 
>       `cd ~/transcriptomics-supbiotech-2023/`
> 
>       `mkdir -p alignment/compressed-results`
> 
>       `samtools sort -n -O BAM alignment/results/T0.sam -o alignment/compressed-results/T0_tmp1.bam`
> 
> - Keep only pairs with both reads mapped:
> 
>       `cd ~/Desktop/transcriptomics-supbiotech-2023/alignment/compressed-results`
> 
>       `samtools fixmate --output-fmt bam T0_tmp1.bam T0_tmp2.bam`
> 
> - Filter reads with low mapping quality:
> 
>       'samtools view --output-fmt bam -f 2 -q 30 -1 -b T0_tmp2.bam -o T0_tmp1.bam'
> 
> - Sort the resulting reads:
> 
>       `samtools sort --output-fmt bam -l 9 T0_tmp1.bam -o T0_sorted.bam`
> 
> - Index the reads
> 
>       `samtools index T0_sorted.bam`
> 
> - Remove the sam files and the temporary files:
> 
>      `rm ../results/T0.sam T0_tmp1.bam T0_tmp2.bam`. 

> The `rm` function allows permanently removing files through the command-line

Repeat for T10.

### Generating the transcription tracks

After the alignment of the reads, it is important to compute the coverage of the genome at each sequence or gene of interest. 

To do that, we use the `bamCoverage` function from `deeptools`.

> 'cd /transcriptomics-supbiotech-2023/alignment/' 
> 
> 'mkdir -p results/tracks'
>
> `bamCoverage --bam compressed-results/T0_sorted.bam --outFileName results/tracks/T0_unstranded.bw --binSize 1 --normalizeUsing CPM --extendReads --ignoreDuplicates`
bamCoverage --bam results/alignment/T0_sorted.bam --outFileName results/tracks/T0_forward.bw --binSize 1 --normalizeUsing CPM --extendReads --filterRNAstrand forward --ignoreDuplicates
bamCoverage --bam results/alignment/T0_sorted.bam --outFileName results/tracks/T0_reverse.bw --binSize 1 --normalizeUsing CPM --extendReads --filterRNAstrand reverse --ignoreDuplicates