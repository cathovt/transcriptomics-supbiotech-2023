## Trimming the adaptors

In the beginning of the reads, there are some **adapter sequences** that need to be removed to avoid mapping issues as they were artificially added.

**_Bacillus_ infection, _T0_**

Read 1: ![Снимок экрана 2023-04-10 в 13.23.28.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar
%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_8ez1hz%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.23.28.png)
Read 2:
![Снимок экрана 2023-04-10 в 13.24.22.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_TD8sRk%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.24.22.png)

_Remove 8 bp:_ 
> cutadapt -u 8 -U 8 -q 33,33 --trim-n -m 20 -o cutadapt/results/Bacillus_infection_T0_R1_trimmed.fq.gz -p cutadapt/results/Bacillus_infection_T0_R2_trimmed.fq.gz fastq/Bacillus_infection_T0_R1.fq.gz fastq/Bacillus_infection_T0_R2.fq.gz

_Check in the command-line:_
> cd cutadapt/results 
> 
> zcat < Bacillus_infection_T0_R1_trimmed.fq.gz | head
> 
> zcat < Bacillus_infection_T0_R2_trimmed.fq.gz | head

**_Bacillus_ infection, _T10_**
Read 1:
![Снимок экрана 2023-04-10 в 13.33.11.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_j74SiL%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.33.11.png)
Read 2:
![Снимок экрана 2023-04-10 в 13.33.37.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_nhBt28%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.33.37.png)

_Remove 8 bp:_
> cutadapt -u 8 -U 8 -q 33,33 --trim-n -m 20 -o cutadapt/results/Bacillus_infection_T10_R1_trimmed.fq.gz -p cutadapt/results/Bacillus_infection_T10_R2_trimmed.fq.gz fastq/Bacillus_infection_T10_R1.fq.gz fastq/Bacillus_infection_T10_R2.fq.gz

_Check in the command-line:_
> cd cutadapt/results 
> 
> zcat < Bacillus_infection_T10_R1_trimmed.fq.gz | head
> 
> zcat < Bacillus_infection_T10_R2_trimmed.fq.gz | head


## Verification
### FastQC

> cd cutadapt   
> fastqc -o results/fastqc results/cutadapt/*fq.gz

**Trimmed _Bacillus_ infection, _T0_**

Read 1:
![Снимок экрана 2023-04-10 в 13.45.32.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_ndY1ON%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.45.32.png)
Read 2:
![Снимок экрана 2023-04-10 в 13.46.03.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_8cJlD8%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.46.03.png)

**Trimmed _Bacillus_ infection, _T10_**

Read 1: 
![Снимок экрана 2023-04-10 в 13.49.26.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_2eCN4N%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.49.26.png)
Read 2:
![Снимок экрана 2023-04-10 в 13.48.43.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_70UWjs%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2013.48.43.png)

### FastQ-Screen

> cd transcriptomics-supbiotech-2023/fastq-screen/FastQ-Screen-0.15.3
> 
> fastq_screen --outdir ../../cutadapt/results/fastq_screen --conf fastq_screen.conf --force ../../cutadapt/results/cutadapt/*fq.gz

**Untrimmed _Bacillus_ infection, _T0_, R1**:
![Снимок экрана 2023-04-10 в 14.00.21.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_ICR9S3%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2014.00.21.png)
![Снимок экрана 2023-04-10 в 14.02.18.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_Ov4hhE%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2014.02.18.png)

Compare with **Trimmed _Bacillus_ infection, _T0_, R1**:
![Снимок экрана 2023-04-10 в 14.06.50.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_b0rrTb%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2014.06.50.png)
![Снимок экрана 2023-04-10 в 14.07.18.png](..%2F..%2F..%2F..%2F..%2F..%2F..%2F..%2Fvar%2Ffolders%2F44%2F_1_t7gpx07l0r0l8gv3b29f80000gn%2FT%2FTemporaryItems%2FNSIRD_screencaptureui_UOP9Y7%2F%D0%A1%D0%BD%D0%B8%D0%BC%D0%BE%D0%BA%20%D1%8D%D0%BA%D1%80%D0%B0%D0%BD%D0%B0%202023-04-10%20%D0%B2%2014.07.18.png)

The fastq-screen analysis results for other reads and time-points are available in: test/results/trimming-adaptors/fastq-screem _**(trimmed)**_ and test/results/fastq-screen **_(untrimmed)_**

## Conclusion: 
We have removed the bias at the beginning of the reads and selected the reads with the highest quality. The length is not uniform now as some reads have been more trimmed than others. Finally, we have resolved the mapping issue from the reverse reads (namely, large fraction of Hit_No_Genomes).