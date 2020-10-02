---
Title: FASTQ to Variant Call (Part 2)
Status: draft
Date: 2020-10-02 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics, genomic variation, entrez-direct, EDirect
---

obtain raw reads from human genomes samples deposited in 

## The FASTQ format

The `FASTQ` file format is how we store the output of whole genome or transcriptomic sequencing (sequences of nucleotides). It inherits its name from the `FASTA` format that stores and the word `Qualities`, because a `FASTQ` file not only contains the nucleotide sequence, but also contains the quality of the sequencing procedure.

The qualities are represented by Phred scores (`Q`), which is used to calculate the probability of a nucleotide being incorrectly identified during sequencing using a formula (I will not go into details here). So, for example, if we check a `FASTQ` file and found a nucleotide with `Q = 30`, it means that there is a probability of 1 in 1000 that it was incorrectly assigned during sequencing -- in other words an accuracy of 99.9%. Therefore, `Q` values around 30 and above are generally seem as very good quality.

### The reason `FASTQ` files contain information about quality

Because during use of these kind of files, it is important that we have confidence on the sequence assignment. During processing in Bioinformatics analysis pipelines, we can remove low-quality nucleotides to ensure that he have the "cleanest" information possible.

### Obtaining FASTQ files

We obtain `FASTQ` after sequencing of genomic samples in platforms such as [Illumina](https://www.illumina.com). Normally, researchers deposit raw `FASTQ` files on public databases to share their discoveries with other scientists. Thus, I will demonstrate how we download publicly-available `FASTQ` files using the recently created `miniconda` environment.

## A Warning

`FASTQ` files, especially from human samples, have very big sizes, in the gigabytes range. Since they are so big, considerable computing power and storage are needed to process these kind of files. Common domestic desktop PCs or laptops are inadequate, because they lack sufficient RAM and CPU processing power to perform smoothly during most Bioinformatic uses.

However, do not despair! I will demonstrate the script using especial `FASTQ` (and other) files that have sufficiently small size so your PC can process it on a timely manner. The logic is: if the script is working in these small files, it should work with the "real deal".
