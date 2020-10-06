---
Title: FASTQ to Annotation (Part 2)
Status: draft
Date: 2020-10-02 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics, genomic variation, entrez-direct, EDirect
---

## Introduction

*In a [previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis), I showed how to configure an Ubuntu system to install Bioinformatics programs.*

*Now, using the environment I created, I will demonstrate a bash script, `FastQ_to_Annotation.sh` that takes next generation sequencing (NGS) raw reads from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.*

*In the [first part](https://antoniocampos13.github.io/fastq-to-annotation-part-1), I showed how to search for NGS projects deposited in [National Center for Biotechnology Information (NCBI) databases](https://www.ncbi.nlm.nih.gov/) from which I can download sequencing reads later to use with the script.*

Here in the second part, I will show how to retrieve raw genome sequencing reads in the form of `FASTQ` files, which are deposited in [SRA](https://www.ncbi.nlm.nih.gov/sra).

But first, let's review what `FASTQ` files are.

## What is the the FASTQ format

The `FASTQ` file format is how we store the output of whole genome or transcriptomic sequencing (sequences of nucleotides). It inherits its name from the `FASTA` format that stores and the word `Qualities`, because a `FASTQ` file not only contains the nucleotide sequence, but also contains the quality of the sequencing procedure.

The qualities are represented by Phred scores (`Q`), which is used to calculate the probability of a nucleotide being incorrectly identified during sequencing using a formula (I will not go into details here). So, for example, if we check a `FASTQ` file and found a nucleotide with `Q = 30`, it means that there is a probability of 1 in 1000 that it was incorrectly assigned during sequencing -- in other words an accuracy of 99.9%. Therefore, `Q` values around 30 and above are generally seem as very good quality.

### The reason `FASTQ` files contain information about quality

Because during use of these kind of files, it is important that we have confidence on the sequence assignment. During processing in Bioinformatics analysis pipelines, we can remove low-quality nucleotides to ensure that he have the "cleanest" information possible.

### Obtaining FASTQ files

We obtain `FASTQ` after sequencing of genomic samples in platforms such as [Illumina](https://www.illumina.com), which practically dominates the NGS market nowadays. Check the fundamentals of Illumina's NGS platform [here](https://www.illumina.com/science/technology/next-generation-sequencing/beginners.html). Normally, researchers deposit raw `FASTQ` files on public databases to share their discoveries with other scientists. This is why I took note of the `BioProject` accession ID during the demonstration of Part 1. With this ID, I can retrieve sequencing reads associated with the project.

### A Warning

`FASTQ` files, especially from human samples, have very big sizes, in the gigabytes range. Therefore, considerable computing power and storage are needed to process these kind of files.

## Retrieving reads from a BioProject

**Activate** the environment, if needed, and connect to `SRA` database via `EDirect esearch` command using the `PRJNA436005` as keyword for query. Then, we pipe the results to the `efetch` command. With the `-format` flag, it will format the results into the `runinfo` format (more on that later). Finally, will save it into the `PRJNA436005_runinfo.csv` file. You can choose other name if you wish.

```bash
# Continuing into the folder I created in the previous part
cd demo
conda activate bioenv

esearch -db sra -query 'PRJNA436005' | efetch -format runinfo > PRJNA436005_runinfo.csv
```

The `runinfo` format displays metadata of read sets. Reads are inferred sequences of base pairs corresponding to DNA fragments produced during procedures for NGS. The collection of DNA fragments from a given sample is called a **library**, which are sequenced to produce the set of **reads**.

Checking the `CSV` file, I can see that there are seven read sets, each displayed on a row, and are identified by the `SRR` prefix followed by some numbers. With this ID is possible to retrieve `FASTQ` files for each read set. Now I check the `LibraryLayout` column to see that they are all **PAIRED** reads, meaning that the researchers sequenced both ends of a fragment. Thus, each read set will produce two `FASTQ`files, containing the sequences and qualities from all reads obtained from the library of the original sample. This is important to check because I will have to input this information in the next commands.

Other interesting columns that I like to check are:

* `spots`, which are the number of physical locations in the sequencing flowcells where the sequencing adapters are fixed. A spot contains several nucleotide bases from several, possibly millions, of reads;
* `avgLength`, which as the name implies, is the average length, in nucleotides, of reads in the set;
* `size_MB`, the size in megabytes of the read set;
* `LibrarySource`, which indicates if the sample source is GENOMIC, TRANSCRIPTOMIC and so on;
* `Platform`, the vendor of NGS procedure;
* `Model`, the model of the `Platform`;
* `Sex`, `Disease` and `Tumor`: descriptors of sample phenotype.

For now, I will use only the first read, which has the `SRR6784104` ID. Finally, let's download the read set with the `EDirect fastq-dump` command and split it into two files (`--split-files` flag), one with reads from each end of DNA fragments in the original library, and compress them with `--gzip`:

```bash
fastq-dump --split-files SRR6784104 --gzip
```

After a moment, two `fastq.gz` files will be downloaded to the current working directory and are ready to be used as the input for the `FastQ_to_VariantCall.sh`.

## Conclusion (Part 2)

In this part I showed how to:

* obtain and inspect metadata from projects via their `BioProjects` accession ID;
* download read sets via their `SRR` accession ID via `fastq-dump`.

## References

[Setting Up Your Unix Computer for Bioinformatics Analysis](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis)

[FASTQ to Annotation (Part 1)](https://antoniocampos13.github.io/fastq-to-annotation-part-1)

[National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)

[Home - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra)

[Illumina | Sequencing and array-based solutions for genetic research](https://www.illumina.com)

[Next-Generation Sequencing for Beginners | NGS basics for researchers](https://www.illumina.com/science/technology/next-generation-sequencing/beginners.html)
