## Introduction

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: rpoetreads from small-RNA sequencing contain the 3’ sequencing adapter because the read is longer than the molecule that is sequenced. Amplicon reads start with a primer sequence. Poly-A tails are useful for pulling out RNA from your sample, but often you don’t want them to be in your reads.
Cutadapt helps with these trimming tasks by finding the adapter or primer sequences in an error-tolerant way. It can also modify and filter single-end and paired-end reads in various ways. Adapter sequences can contain IUPAC wildcard characters. Cutadapt can also demultiplex your reads.
Cutadapt is available under the terms of the MIT license.
Cutadapt development was started at TU Dortmund University in the group of Prof. Dr. Sven Rahmann. It is currently being developed within NBIS (National Bioinformatics Infrastructure Sweden).
If you use Cutadapt, please cite DOI:10.14806/ej.17.1.200 .
## Installation
>pip install cutadapt

[PyPi page](https://pypi.org/project/cutadapt/)

