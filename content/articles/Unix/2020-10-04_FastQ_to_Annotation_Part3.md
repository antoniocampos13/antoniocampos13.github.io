---
Title: FASTQ to Annotation (Part 3)
Status: published
Date: 2020-10-05 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics, genomic variation, entrez-direct, EDirect
---

## Introduction

*In a [previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html), I showed how to configure an Ubuntu system to install Bioinformatics programs.*

*Now, using the environment I created, I will demonstrate a bash script, `FastQ_to_Annotation.sh` that takes next generation sequencing (NGS) raw reads from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.*

*In the [first part](https://antoniocampos13.github.io/fastq-to-annotation-part-1), I showed how to search for NGS projects deposited in [National Center for Biotechnology Information (NCBI) databases](https://www.ncbi.nlm.nih.gov/) from which I can download sequencing reads later to use with the script.*

*In the [second part](https://antoniocampos13.github.io/fastq-to-annotation-part-2), I showed how to retrieve raw genome sequencing reads in the form of `FASTQ` files, which are deposited in [SRA](https://www.ncbi.nlm.nih.gov/sra).*

Here in the third part, I make the final preparations for the `FastQ_to_Annotation.sh` script demonstration using the `FASTQ` files obtained in the previous part.

## Final preparations

### Installing local cache of Ensembl Variant Effect Predictor (VEP)

The [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) is the core tool used by the script for the annotation of the effects of any variants present in the sample. It may be used online, but Ensembl strongly recommends users to download and install a local cache of all data deposited in the tool to avoid server overload. I open a terminal, activate the `bioenv` miniconda environment and execute the commands below:

**WARNING: Several gigabytes of data will be downloaded from the internet and installed on your computer. Be sure that you have plenty or unlimited data allowances from your ISP and sufficient free space on your hard drive before continuing. It will take a while (several minutes to hours) until all the needed processes finish.**

```bash
vep_install -a cf -s homo_sapiens_refseq -y GRCh38 -c . –CONVERT
```

Unfortunately, the command above is prone to network and other esoteric errors. If you have any problem, you can try an alternative manner. See below:

```bash
# Alternative: manually download the compressed cache
wget ftp://ftp.ensembl.org/pub/release-101/variation/indexed_vep_cache/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -P $HOME/.vep

# Uncompress the cache
tar -zxf $HOME/.vep/homo_sapiens_refseq_vep_101_GRCh38.tar.gz -C $HOME/.vep
```

After some time (several minutes to some hours depending on the network and the capabilities of your computer) the VEP cache will be downloaded and installed into a hidden folder in your home folder (`$HOME/.vep`). Therefore, notice that it is independent of miniconda environments. Thus, it is expected that once installed, this cache will work with any miniconda environment on your computer. You may backup the `.vep/` folder to avoid downloading the whole thing again (but consider to download newer versions of the cache as they become available though).

### Obtaining human genome reference files

The annotation process require a collection of reference files. These files will assist us to generate a "list" of genetic variants, alongside their possible effects on the phenotype (the **annotation**).

These files are:

1. A [`FASTA`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp) file (extensions `.fasta`, `.fa` or `.fna`). It must contain the complete nucleotide sequence of the human genome. We will compare our `FASTQ` files against it;
2. A [`FASTA index`](http://www.htslib.org/doc/faidx.html) (`.fai`). It stores genomic regions as coordinates. We will use it to generate a Browser Extensible Data (`.bed`) file (see below);
3. A [Browser Extensible Data (`.bed`)](https://en.wikipedia.org/wiki/BED_(file_format)) file. It stores genomic regions as coordinates, indicating the start and end of chromosomes. It is most useful when its information is chromosome-ordered and position-sorted;
4. Alternatively, the `FastQ_to_Annotation.sh` script can accept a [General Feature Format (`.gff`)](https://www.ensembl.org/info/website/upload/gff.html) instead of the `.bed` file. It is used for describing genes and other features of DNA, RNA and protein sequences;
5. Burrows-Wheelers Aligner index files (`.amb`, `.ann`, `.bwt`, `.pac` and `.sa`). The [`bwa` program](http://bio-bwa.sourceforge.net/) is a short read alignment tool. In other words, it identifies the location of the reads inside the `FASTQ` files. The `FastQ_to_Annotation.sh` script uses `bwa` at the second step of the pipeline. It works by efficiently using this collection of five files as a index.

*Why we need these files?* Briefly, They serve to map the genomic location of any variant we identify in our samples, as well as the genetic mutation that occurred there, which allows us to predict the possible effect(s) over the phenotype in question (in our case, MAPKi-resistance in melanoma samples). If we compare the genetic variation profile of MAPKi-susceptible samples with MAPKi-resistant samples, we could identify genetic variants associated with the resistances, and perhaps point to new directions of prognosis and new treatments.

I will now show how to obtain all these files. Remember to **activate** the `miniconda` that you created before if necessary.

**WARNING: Several gigabytes of data will be downloaded from the internet and installed on your computer. Be sure that you have plenty or unlimited data allowances from your ISP and sufficient free space on your hard drive before continuing. It will take a while (several minutes to hours) until all the needed processes finish.**

#### 1. Human genome FASTA

```bash
# First, create a subfolder into the demo folder to better organize our reference files
cd demo
mkdir refs
cd refs

# Download GRCh38 major release without ALT contigs and with decoy genomes (EBV and hs38d1 contig) from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
```

#### 2. Human genome FASTA index

```bash
# Download from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai
```

#### 3. Human genome BED file

```bash
# Produce sorted BED file from reference genome index file obtained above
awk '{print $1 "\t0\t" $2}' GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai | sort -k1,1V -k2,2n > GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.bed
```

#### 4. Human genome GFF file (optional alternative to BED file)

```bash
# Download from NCBI's FTP server
curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
```

#### 5. bwa index files

```bash
bwa index GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
```

When you complete all of the steps in Part 1, Part 2 and in this part, your `demo` folder should have the files showed below

![demo folder contents until now]({static}/images/demo_folder.PNG)

Now, go to the folder [`FastQ_to_Annotation` folder in my portfolio](https://github.com/antoniocampos13/portfolio/tree/master/Unix/2020-10-01_Fastq%20to%20Annotation), take heed of the GPL License and Copyright Notice, download and copy the `FastQ_to_Annotation.sh` script into your `demo` folder.

Thus, the only mandatory files are the `FastQ_to_Annotation.sh`, the `FASTQ` pair and the ones in the `refs` folder. If you are missing any other file, do not worry.

## GPL License and Copyright Notice

The `FastQ_to_Annotation.sh` script is a modified version from [Dr. Kevin Blighe's original scripts](https://github.com/kevinblighe/ClinicalGradeDNAseq). Both works are licensed under [GNU General Public License v3.0](https://www.gnu.org/licenses/licenses.en.html).

## Conclusion of Part 3

In this part I showed how to:

* Set up VEP local cache;
* Obtain human genome reference files;
* Obtain auxiliary files needed for short read alignment.

Finally, everything is in place for the `FastQ_to_Annotation.sh` script.

*[Go to FASTQ to Annotation (Part 4)](https://antoniocampos13.github.io/fastq-to-annotation-part-4)*

*[Go back to FASTQ to Annotation (Part 1)](https://antoniocampos13.github.io/fastq-to-annotation-part-1)*

*[Go back to FASTQ to Annotation (Part 2)](https://antoniocampos13.github.io/fastq-to-annotation-part-2)*

## References

[National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/)

[Home - SRA - NCBI](https://www.ncbi.nlm.nih.gov/sra)

[Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)

[BLAST Topics | Query Input and database selection](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp)

[faidx(5) manual page](http://www.htslib.org/doc/faidx.html)

[BED (file format) - Wikipedia](https://en.wikipedia.org/wiki/BED_(file_format))

[GFF/GTF File Format](https://www.ensembl.org/info/website/upload/gff.html)

[Burrows-Wheeler Aligner](http://bio-bwa.sourceforge.net/)

[kevinblighe/ClinicalGradeDNAseq](https://github.com/kevinblighe/ClinicalGradeDNAseq)

[GNU General Public License](https://www.gnu.org/licenses/licenses.en.html)

[NCBI's FTP Server | GRCh38 Major release sequences for alignment pipelines](http://ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/)

[Ensembl's FTP Server | Homo sapiens DNA sequences release 101](http://ftp://ftp.ensembl.org/pub/release-101/fasta/homo_sapiens/dna/)
