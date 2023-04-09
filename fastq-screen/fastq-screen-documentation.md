# Introduction
FastQ Screen is a simple application which allows to search a large sequence dataset against a panel of different genomes to determine from where the sequences in the data originate. It was built as a QC check for sequencing pipelines but may also be useful to characterise metagenomic samples. When running a sequencing pipeline it is useful to know that sequencing runs contain the types of sequence they’re supposed to. Your search libraries might contain the genomes of all of the organisms you work on, along with PhiX, Vectors or other contaminants commonly seen in sequencing experiments.

Although the program wasn’t built with any particular technology in mind it is probably only really suitable for processing short reads due to the use of either Bowtie, Bowtie2 or BWA as the searching application.

The program generates both text and graphical output to inform you what proportion of your library was able to map, either uniquely or to more than one location, against each of your specified reference genomes. The user should therefore be able to identify a clean sequencing experiment in which the overwhelming majority of reads are probably derived from a single genomic origin.
# FastQ Screen online tutorials
[Introduction to FastQ Screen](https://www.youtube.com/watch?v=8IsGdikLhaE)

[Downloading, configuring and running FastQ Screen](https://www.youtube.com/watch?v=WqiKPRxHzNU)

[Interpreting FastQ Screen results](https://www.youtube.com/watch?v=x32k84HHqjQ)

[Filtering FASTQ files](https://www.youtube.com/watch?v=eJcAv-Dt57I)

# Download
FastQ Screen may be obtained from the [GitHub download page](https://github.com/StevenWingett/FastQ-Screen/releases/tag/v0.15.3).
# Installation
Before running FastQ Screen there are a few prerequisites that will need to be installed:

- A sequence aligner. FastQ Screen is compatible with **Bowtie**, **Bowtie2** or **BWA**. It’s easier if you put the chosen aligner in your path, but if not you can configure its location in the config file.

- We recommend running FastQ Screen in a Linux system, on which the programming language Perl should already be installed.

- **GD::Graph**. FastQ Screen uses the GD::Graph module to draw PNG format graphs summarising the mapping results. FastQ Screen will still produce both text and HTML format summaries of the results if GD::Graph is not installed.

### Configuration
Modify fastq_screen.conf file to specify the path to aligner (Bowtie / Bowtie 2 / BWA) and add genome databases.

Example of a configuration:
Used configuration: fastq_screen.conf

_Example of an execution command:_
> FastQ-Screen-0.15.3/fastq_screen --outdir results --conf FastQ-Screen-0.15.3/fastq_screen.conf --force ../fastq/Bacillus_infection_T0_R1.fq.gz ../fastq/Bacillus_infection_T0_R2.fq.gz .//fastq/Bacillus_infection_T10_R1.fq.gz .//fastq/Bacillus_infection_T10_R2.fq.gz
