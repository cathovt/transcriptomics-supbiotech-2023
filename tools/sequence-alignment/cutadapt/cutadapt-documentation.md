## Introduction

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Cleaning your data in this way is often required: reads from small-RNA sequencing contain the 3’ sequencing adapter because the read is longer than the molecule that is sequenced. Amplicon reads start with a primer sequence. Poly-A tails are useful for pulling out RNA from your sample, but often you don’t want them to be in your reads.

Cutadapt helps with these trimming tasks by finding the adapter or primer sequences in an error-tolerant way. It can also modify and filter single-end and paired-end reads in various ways. Adapter sequences can contain IUPAC wildcard characters. Cutadapt can also demultiplex your reads.

Cutadapt is available under the terms of the MIT license.

Cutadapt development was started at TU Dortmund University in the group of Prof. Dr. Sven Rahmann. It is currently being developed within NBIS (National Bioinformatics Infrastructure Sweden).

## Installation
>pip install cutadapt

[PyPi page](https://pypi.org/project/cutadapt/)

See full documentation [here](https://cutadapt.readthedocs.io/).

To see all command-line options:
> cutadapt --help"

## Options:

>  _-h, --help_:         show this help message and exit
>
>  _--version_:           show version number and exit
>
>  _--debug_:             print debug log. Use twice to also print DP matrices
> 
> _-j CORES, --cores CORES_: number of CPU cores to use. Use 0 to auto-detect. Default: 1

## Finding adapters:
  _Parameters **-a, -g, -b** specify adapters to be removed from each read (or from R1 if data is paired-
  end. If specified multiple times, only the best matching adapter is trimmed (but see the --times
  option). Use notation 'file:FILE' to read adapter sequences from a FASTA file._

>  _-a ADAPTER, --adapter ADAPTER_
>
> Sequence of an adapter ligated to the 3' end (paired data: of the first read). 

The adapter and subsequent bases are trimmed. If a '$' character is appended
                        ('anchoring'), the adapter is only found if it is a suffix of the read. 

>_-g ADAPTER, --front ADAPTER_
>
> Sequence of an adapter ligated to the 5' end (paired data: of the first read).

The adapter and any preceding bases are trimmed. Partial matches at the 5' end
                        are allowed. If a '^' character is prepended ('anchoring'), the adapter is only
                        found if it is a prefix of the read. 

>_-b ADAPTER, --anywhere ADAPTER_
>
> Sequence of an adapter that may be ligated to the 5' or 3' end (paired data: of
                        the first read). Both types of matches as described under -a and -g are allowed.

If the first base of the read is part of the match, the behavior is as with -g,
                        otherwise as with -a. This option is mostly for rescuing failed library
                        preparations - do not use if you know which end your adapter was ligated to!

> _-e E, --error-rate E, --errors E_

Maximum allowed error rate (if 0 <= E < 1), or absolute number of errors for
                        full-length adapter match (if E is an integer >= 1). Error rate = no. of errors
                        divided by length of matching region. Default: 0.1 (10%)

> _--no-indels_          
>
> Allow only mismatches in alignments. Default: allow both mismatches and indels 

> _-n COUNT, --times COUNT_
>
> Remove up to COUNT adapters from each read. Default: 1

> _-O MINLENGTH, --overlap MINLENGTH_
>
> Require MINLENGTH overlap between read and adapter for an adapter to be found. Default: 3

> _--match-read-wildcards_
>
> Interpret IUPAC wildcards in reads. Default: False

> _-N, --no-match-adapter-wildcards_
>
> Do not interpret IUPAC wildcards in adapters.

> _--action {trim,retain,mask,lowercase,none}_
> What to do if a match was found. 
> 
> _trim:_ trim adapter and up- or downstream sequence; 
> 
> _retain_: trim, but retain adapter;
> 
> _mask_: replace with 'N' characters; 
> 
> _lowercase_: convert to lowercase; 
> 
> _none_: leave unchanged. 
> 
> **Default**: _trim_

> --rc, --revcomp       
Check both the read and its reverse complement for adapter matches. If match is
                        on reverse-complemented version, output that one. Default: check only read

## Additional read modifications:

> _-u LEN, --cut LEN_
> 
> Remove LEN bases from each read (or R1 if paired; use -U option for R2). If LEN
                        is positive, remove bases from the beginning. If LEN is negative, remove bases
                        from the end. Can be used twice if LENs have different signs. Applied *before*
                        adapter trimming.

> _--nextseq-trim 3'CUTOFF_
> 
> NextSeq-specific quality trimming (each read). Trims also dark cycles appearing
                        as high-quality G bases.

> _-q [5'CUTOFF,]3'CUTOFF, --quality-cutoff [5'CUTOFF,]3'CUTOFF_
>
> Trim low-quality bases from 5' and/or 3' ends of each read before adapter
                        removal. Applied to both reads if data is paired. If one value is given, only
                        the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is
                        trimmed with the first cutoff, the 3' end with the second.
> _--quality-base N_     
>
> Assume that quality values in FASTQ are encoded as ascii(quality + N). This
                        needs to be set to 64 for some old Illumina FASTQ files. 
> 
> **Default:** 33

> _--length LENGTH, -l LENGTH_
>
> Shorten reads to LENGTH. Positive values remove bases at the end while negative
                        ones remove bases at the beginning. This and the following modifications are
                        applied after adapter trimming.
> _--trim-n_              
>
> Trim N's on ends of reads. 

> _length-tag TAG_      
>
>Search for TAG followed by a decimal number in the description field of the
                        read. Replace the decimal number with the correct length of the trimmed read.
> 
> _For example_, use --length-tag 'length=' to correct fields like 'length=123'.

> --strip-suffix STRIP_SUFFIX
> 
> Remove this suffix from read names if present. Can be given multiple times.

> -_x PREFIX, --prefix PREFIX_
>
> Add this prefix to read names. Use {name} to insert the name of the matching
                        adapter.

> _-y SUFFIX, --suffix SUFFIX_
>
> Add this suffix to read names; can also include {name}

> _--rename TEMPLATE_     
>
> Rename reads using TEMPLATE containing variables such as {id}, {adapter_name}
                        etc. (see documentation)
> _--zero-cap, -z_       
>
> Change negative quality values to zero.

### Filtering of processed reads:
  _Filters are applied after above read modifications. Paired-end reads are always discarded pairwise
  (see also --pair-filter)._

>  _-m LEN[:LEN2], --minimum-length LEN[:LEN2]_     
> Discard reads shorter than LEN. **Default**: 0
> 
> -_M LEN[:LEN2], --maximum-length LEN[:LEN2]_  
> Discard reads longer than LEN. **Default**: no limit
> 
> _--max-n COUNT_         
> Discard reads with more than COUNT 'N' bases.
> If COUNT is a number between 0 and 1, it is interpreted as a fraction of the read length.
>
> _--max-expected-errors ERRORS, --max-ee ERRORS_       
> Discard reads whose expected number of errors (computed from quality values)
                        exceeds ERRORS.
>
> _--discard-trimmed, --discard_      
> Discard reads that contain an adapter. Use also -O to avoid discarding too many
                        randomly matching reads.
>
> _--discard-untrimmed, --trimmed-only_       
Discard reads that do not contain an adapter.
> 
> _--discard-casava_      
> Discard reads that did not pass CASAVA filtering (header has :Y:).

### Output:
>  _--quiet_        
> Print only error messages.

> _--report {full,minimal}_        
> Which type of report to print: 'full' or 'minimal'. **Default**: full 

> _--json FILE_     
> Dump report in JSON format to FILE 

> _-o FILE, --output FILE_      
> Write trimmed reads to FILE. FASTQ or FASTA format is chosen depending on input.
                        Summary report is sent to standard output. Use '{name}' for demultiplexing (see
                        docs). **Default:** write to standard output

> _--fasta_       
> Output FASTA to standard output even on FASTQ input.

> _-Z_        
> Use compression level 1 for gzipped output files (faster, but uses more space)

> _--info-file FILE_             
> Write information about each read and its adapter matches into FILE. See the
                        documentation for the file format. 

> _-r FILE, --rest-file FILE_       
> When the adapter matches in the middle of a read, write the rest (after the
                        adapter) to FILE.

> _--wildcard-file FILE_        
> When the adapter has N wildcard bases, write adapter bases matching wildcard
                        positions to FILE. (Inaccurate with indels.)

> _--too-short-output FILE_     
Write reads that are too short (according to length specified by -m) to FILE.
                        Default: discard reads 

> _--too-long-output FILE_      
> Write reads that are too long (according to length specified by -M) to FILE.
                        Default: discard reads 

> _--untrimmed-output FILE_     
> Write reads that do not contain any adapter to FILE. Default: output to same
                        file as trimmed reads

### Paired-end options:
  _The -A/-G/-B/-U/-Q options work like their lowercase counterparts, but are applied to R2 (second
  read in pair)_

>  _-A ADAPTER_     
> 3' adapter to be removed from R2
>   
> _-G ADAPTER_        
> 5' adapter to be removed from R2
>
> _-B ADAPTER_       
> 5'/3 adapter to be removed from R2
> 
> _-U LENGTH_       
> Remove LENGTH bases from R2
>
> _-Q [5'CUTOFF,]3'CUTOFF_      
> Quality-trimming cutoff for R2. **Default**: same as for R1
>
> _-p FILE, --paired-output FILE_       
> Write R2 to FILE.
>
> _--pair-adapters_     
> Treat adapters given with -a/-A etc. as pairs. Either both or none are removed
                        from each read pair.
> 
> _--pair-filter {any,both,first}_      
> Which of the reads in a paired-end read have to match the filtering criterion in
                        order for the pair to be filtered. **Default**: any
> 
> _--interleaved_       
> Read and/or write interleaved paired-end reads.
>
> _--untrimmed-paired-output FILE_      
 Write second read in a pair to this FILE when no adapter was found. Use with
_--untrimmed-output_. **Default**: output to same file as trimmed reads
>
> _--too-short-paired-output FILE_         
> Write second read in a pair to this file if pair is too short.
> 
> _--too-long-paired-output FILE_     
> Write second read in a pair to this file if pair is too long.
