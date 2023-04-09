# Overview of FastQC Results

- Longer forward reads (51bp) than reverse reads (35bp)

**Forward read (51bp):**

<img width="878" alt="Снимок экрана 2023-04-09 в 13 06 02" src="https://user-images.githubusercontent.com/56762339/230766686-bebc203e-59e8-4c96-9209-e6c10c652dab.png">

**Reverse read (35bp):**

<img width="890" alt="Снимок экрана 2023-04-09 в 13 12 44" src="https://user-images.githubusercontent.com/56762339/230766985-796f7d0d-635d-4559-9041-9b4470159fb8.png">


- No issue of sequence quality

**Forward read (51bp):**

<img width="846" alt="Снимок экрана 2023-04-09 в 13 06 54" src="https://user-images.githubusercontent.com/56762339/230766723-8064670b-ae70-4c2f-b724-3cf900b7fe21.png">

**Reverse read (35bp):**

<img width="863" alt="Снимок экрана 2023-04-09 в 13 13 01" src="https://user-images.githubusercontent.com/56762339/230766995-37ba25b1-2f09-441f-b7a7-91bc5bea3371.png">

- We have an expected biased fragmentation, especially strong in reverse reads which is expected in RNAseq

**Forward read (51bp):**

<img width="885" alt="Снимок экрана 2023-04-09 в 13 07 49" src="https://user-images.githubusercontent.com/56762339/230766759-fd802490-4f6b-43ea-aba9-064e0a9557b1.png">

**Reverse read (35bp):**

<img width="869" alt="Снимок экрана 2023-04-09 в 13 13 38" src="https://user-images.githubusercontent.com/56762339/230767018-c6cb66d0-6293-4630-9f34-418b364b5bdd.png">

- Weird distributions in Reverse GC content with a scale stuff, might be due to the contamination with adapters

**Forward read (51bp):**

<img width="900" alt="Снимок экрана 2023-04-09 в 13 10 26" src="https://user-images.githubusercontent.com/56762339/230766872-4e09c04f-c1c5-43cf-bcc2-b54080e3a7d2.png">

**Reverse read (35bp):**

<img width="878" alt="Снимок экрана 2023-04-09 в 13 14 09" src="https://user-images.githubusercontent.com/56762339/230767043-15c0ce5d-19b1-4f09-b94d-e65ec0a75b0f.png">

- No N sequences, fine read length

**Forward read (51bp):**

<img width="866" alt="Снимок экрана 2023-04-09 в 13 10 56" src="https://user-images.githubusercontent.com/56762339/230766892-2d7561a3-5d18-42fb-ac3f-953cd7efbe9a.png">

**Reverse read (35bp):**

<img width="873" alt="Снимок экрана 2023-04-09 в 13 14 21" src="https://user-images.githubusercontent.com/56762339/230767056-2a92fbd0-3454-413a-8fd9-e12f8607139c.png">

- Small duplication levels, expected for RNA-seq

**Forward read (51bp):**

<img width="933" alt="Снимок экрана 2023-04-09 в 13 11 19" src="https://user-images.githubusercontent.com/56762339/230766918-da950a88-0f59-4e98-bab4-e3cb1c2cb40a.png">

**Reverse read (35bp):**

<img width="845" alt="Снимок экрана 2023-04-09 в 13 14 49" src="https://user-images.githubusercontent.com/56762339/230767076-76502790-3843-4fbe-8f70-c28b207e81dd.png">

- No adapters content

**Forward read (51bp):**

<img width="878" alt="Снимок экрана 2023-04-09 в 13 12 12" src="https://user-images.githubusercontent.com/56762339/230766962-1ac99aed-d313-4422-95ba-f43b5072d017.png">

**Reverse read (35bp):**

<img width="874" alt="Снимок экрана 2023-04-09 в 13 15 40" src="https://user-images.githubusercontent.com/56762339/230767111-fae17a59-dd00-49ed-bc30-babfeb3b5967.png">

- Small overepresented GGGGGGGGGGGG reads in reverse librairies which is due to insufficient DNA output

**Forward read (51bp):**

<img width="799" alt="Снимок экрана 2023-04-09 в 13 11 43" src="https://user-images.githubusercontent.com/56762339/230766947-2ee5a71e-72a4-4ab4-850f-2becd61a0f44.png">

**Reverse read (35bp):**

<img width="626" alt="Снимок экрана 2023-04-09 в 13 15 07" src="https://user-images.githubusercontent.com/56762339/230767088-3a5ccf96-2c84-4394-881e-d182ebfc395a.png">
