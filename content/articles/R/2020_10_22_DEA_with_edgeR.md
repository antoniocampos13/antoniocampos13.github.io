---
Title: Differential Expression Analysis with edgeR in R
Status: draft
Date: 2020-10-22 15:30
Author: Antonio Victor Campos Coelho
Categories: R
Tags: Bioinformatics, gene expression, edgeR
---

## Introduction

In my [previous post](https://antoniocampos13.github.io/data-manipulation-with-r.html#data-manipulation-with-r) I demonstrated how to organize the [CGC prostate cancer data](https://antoniocampos13.github.io/working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1.html) to a format suited to differential expression analysis (DEA).

Nowadays, DEA usually arises from high-throughput sequencing of a collection (library) of RNA molecules expressed by single cells or tissue given their conditions upon collection and RNA extraction.

In terms of statistical analysis, DEA ["means taking the normalized read count data and performing statistical analysis to discover quantitative changes in expression levels between experimental groups"](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene). What are experimental groups? Consider for example, diseased versus healthy cells, treated cells versus non-treated cells (when someone is testing new drugs for example), and so on.

## The edgeR package

There are some statistical packages in R that deal with DEA, such as `edgeR`, `DESeq2` and `limma`. Here I will demonstrate a custom script to perform DEA with [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html).

`edgeR` performs DEA for pre-defined genomic features, which can be genes, transcripts or exons, for example. In the present demonstration, we will quantify transcripts. Remember that genes can produce several transcripts through alternative splicing. The `edgeR` statistical model is based on negative binomial distribution. Prior to statistical analysis, `edgeR` normalizes gene/transcript expression counts via the Trimmed Mean of M-values (TMM) method. See [Dr. Kevin Blighe's comment in a Biostars forum topic](https://www.biostars.org/p/284775/#284893) for a brief discussion of `edgeR` and other DEA packages.

Without further ado, I will show how to set up a R session to run a DEA with `edgeR`, and how to interpret results. As usual, the code presented here is deposited on my [portolio at GitHub]().

## Install and load packages

First, I will install two packages. The first one is [`openxlsx`](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.2), which is a package used to read/write Microsoft Office Excel spreadsheets. I will use it to conveniently save the output of the DEA.

The second is [`BiocManager`](https://www.rdocumentation.org/packages/BiocManager/versions/1.30.10). It is needed to install packages from [Bioconductor](https://www.bioconductor.org/) project. They are Bioinformatics analysis packages that are not on the default R package repository.

```r
# Run only once
install.packages(c("openxlsx", "BiocManager"))
```

Now I install `edgeR` and some more packages from Bioconductor. I will use them to annotate and convert the transcript/gene IDs to a gene symbol. Check their documentation: [`annotate`](https://www.rdocumentation.org/packages/annotate/versions/1.50.0), [`org.Hs.eg.db`](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html), [EnsDb.Hsapiens.v79](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v79.html) and [`ensembldb`](https://www.rdocumentation.org/packages/ensembldb/versions/1.4.7):

```r
# Run only once
BiocManager::install(c("edgeR", "annotate", "org.Hs.eg.db", "EnsDb.Hsapiens.v79","ensembldb"))
```

The `::` is used when we wish to invoke the mentioned package directly (`BiocManager` in this case), without loading it into the memory.

Now, I will load just the `here` package to handle filepaths for now, the rest will be loaded into R later (I am assuming you installed it from last time).

```r
library(here)
```

I will use the `counts` data frame I produced [last time](https://antoniocampos13.github.io/data-manipulation-with-r.html). Since I have saved it to my disk, I load it into the current R session:

```r
load(here("data", "counts.RData"))
```
