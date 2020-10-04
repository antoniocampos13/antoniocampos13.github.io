---
Title: FASTQ to Variant Call (Part 3)
Status: draft
Date: 2020-10-05 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics, genomic variation, entrez-direct, EDirect
---

## Introduction

*In a [previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html), I showed how to configure an Ubuntu system to install Bioinformatics programs.*

*Now, using the environment I created, I will demonstrate a bash script, `FasQ_to_VariantCall.sh` that takes next generation sequencing (NGS) raw reads from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.*

*In the first part, I showed how to search for NGS projects deposited in [National Center for Biotechnology Information (NCBI) databases](https://www.ncbi.nlm.nih.gov/) from which I can download sequencing reads later to use with the script.*

*In the second part, I showed how to retrieve raw genome sequencing reads in the form of `FASTQ` files, which are deposited in [SRA](https://www.ncbi.nlm.nih.gov/sra).*

Here in the third and final part, everything is ready for the `FasQ_to_VariantCall.sh` script demonstration using the `FASTQ` files obtained in the [previous part](https://antoniocampos13.github.io/fastq-to-variant-call-part-2).

## Final preparations

### Installing local cache of Ensembl Variant Effect Predictor (VEP)

The [Ensembl Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html) is the core tool used by the script for the annotation of the effects of any variants present in the sample. It may be used online, but Ensembl recommends to download and install a local cache of all data deposited in the tool to avoid server overload.

**Activate** the miniconda environment if needed and execute the command:

**WARNING** about 100 GB of data will be downloaded from the internet and installed on your computer. Be sure that you have plenty or unlimited data allowances from your ISP and sufficient free space on your hard drive before continuing.

```bash
vep_install -a cf -s homo_sapiens_refseq -y GRCh38 -c . –CONVERT
```

By default, the cache will be installed into a hidden folder in your home folder (`/home/<YOUR_USER_NAME>/.vep`). Thus, it is expected that different miniconda environments on the same system can find this cache. You may backup the `.vep/` folder to avoid downloading the whole thing again (but consider to download newer versions of the cache as they become available though). Remember to restore it to the equivalent location as needed (again: `/home/<YOUR_USER_NAME>/.vep`).

### Ensuring that the script will use the correct miniconda environment

## GPL License and Copyright Notice

The `FasQ_to_VariantCall.sh` script is a modified version from [Dr. Kevin Blighe's original scripts](https://github.com/kevinblighe/ClinicalGradeDNAseq). Both works are licensed under [GNU General Public License v3.0](http://www.gnu.org/licenses).
