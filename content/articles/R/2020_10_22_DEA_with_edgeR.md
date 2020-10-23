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

There are some statistical packages in R that deal with DEA, such as `edgeR`, `DESeq2` and `limma`. Here I will demonstrate a custom script to perform DEA with [`edgeR`](https://bioconductor.org/packages/release/bioc/html/edgeR.html). The demonstration here is on Windows 10, but the same steps can be performed on Unix systems.

`edgeR` performs DEA for pre-defined genomic features, which can be genes, transcripts or exons, for example. In the present demonstration, we will quantify transcripts. Remember that genes can produce several transcripts through alternative splicing. The `edgeR` statistical model is based on negative binomial distribution. Prior to statistical analysis, `edgeR` normalizes gene/transcript expression counts via the Trimmed Mean of M-values (TMM) method. See [Dr. Kevin Blighe's comment in a Biostars forum topic](https://www.biostars.org/p/284775/#284893) for a brief discussion of `edgeR` and other DEA packages.

Without further ado, I will show how to set up a R session to run a DEA with `edgeR`, and how to interpret results. As usual, the code presented here is deposited on my [portfolio at GitHub](https://github.com/antoniocampos13/portfolio/tree/master/R/2020_10_22_DEA_with_edgeR).

## Install and load packages

First, I will install some new packages that I have not talked about. The first one is [`openxlsx`](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.2), which is a package used to read/write Microsoft Office Excel spreadsheets. I will use it to conveniently save the output of the DEA.

The second is [`BiocManager`](https://www.rdocumentation.org/packages/BiocManager/versions/1.30.10). It is needed to install packages from [Bioconductor project](https://www.bioconductor.org/), which hosts Bioinformatics analysis packages that are not on the default R package repository.

```r
# Run only once
install.packages(c("here", "tidyverse", "openxlsx", "BiocManager"))
```

Now I install `edgeR` and some more packages from Bioconductor. I will use them to annotate and convert the transcript/gene IDs to a gene symbol. Check their documentation: [`AnnotationDbi`](https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html), [`annotate`](https://www.rdocumentation.org/packages/annotate/versions/1.50.0), [`org.Hs.eg.db`](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html), [EnsDb.Hsapiens.v79](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v79.html) and [`ensembldb`](https://www.rdocumentation.org/packages/ensembldb/versions/1.4.7):

```r
# Run only once
BiocManager::install(c("edgeR","AnnotationDbi", "annotate", "org.Hs.eg.db", "EnsDb.Hsapiens.v79","ensembldb"))
```

The `::` is used when we wish to invoke the mentioned package directly (`BiocManager` in this case), without loading it into the memory.

Now, I will load just the `here` package to handle file paths for now, the rest will be loaded into R later.

```r
library(here)
```

I will use the `counts` data frame I produced [last time](https://antoniocampos13.github.io/data-manipulation-with-r.html). Since I have saved it to my disk, I load it into the current R session. If you already have the `counts` data frame loaded in the session from the previous demonstration, this step is not necessary.

```r
load(here("data", "counts.RData"))
```

Now I load the custom `edgeR_setup()` function I use to perform DEA with `edgeR`. I wrote the function in a R script with the same name and saved it on my `src` folder:

```r
source(here("src", "edgeR_setup.R"))
```

Check the `edgeR_setup.R` script. First, it loads the packages I installed before:

```r
library(tidyverse)
library(openxlsx)

# Install trough BiocManager
library(edgeR)
library(annotate)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v79)
library(ensembldb)
```

Now, check the function arguments:

```r
edger_setup <- function(name, counts, replicates = TRUE, filter = TRUE, gene_id = c("NCBI", "ENSEMBL", "SYMBOL"), output_path) {

    # ... the function goes here ...
}
```

* `name`: A string. An identifier for the experiment.
* `counts`: The data frame containing the transcript counts.
* `replicates`: A Boolean indicating if the samples are biological replicates. Defaults to `TRUE`.
* `filter`: A Boolean indicating if lowly expressed transcripts should be filter out. Defaults to `TRUE`.
* `gene_id`: A string indicating how transcripts are identified in the data frame. There are three options:
    * `NCBI`: [Entrez Gene Record ids](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/Genes/sample_entrez_gene_record.html)
    * `ENSEMBL`: [ENSEMBL ids (ENS#)](https://www.ebi.ac.uk/training-beta/online/courses/ensembl-browsing-genomes/navigating-ensembl/investigating-a-gene/)
    * `SYMBOL`: [Official HGNC gene symbol](https://www.genenames.org/)
* `output_path`: A path and filename string where the results will be saved in Excel spreadsheet format. Example: `"\some\path\results.xlsx"`.

The use of `edgeR` to analyze datasets with no biological replicates (`replicates = FALSE`) is discouraged. However, I prepared a special dataset of housekeeping genes based on the work by [Eisenberg and Levanon (2013)](https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00089-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0168952513000899%3Fshowall%3Dtrue). I downloaded the `HK_genes.txt` supplementary file [here](https://www.tau.ac.il/~elieis/HKG/) and placed it in the `data` folder at my current work directory. I also wrote an auxiliary script named `hk_genes.R` and placed it in the `src` folder. Check below a representation of my current work directory, where `main_dea_edgeR.R` contains the commands of this demonstration.

```txt
.
├── data
│   ├── counts.RData
│   └── HK_genes.txt
├── output
├── src
│   ├── edgeR_setup.R
│   └── hk_genes.R
└── main_dea_edgeR.R
```

For now, I will use default values and indicate that the transcripts in the data frame are identified by ENSEMBL ids. Before running the function, I assign the output path string to the `out_path` object, which I include in the function call. Note that I gave the name `prostate_cancer` to identify the experiment.

```r
out_path <- here("output", "prostate_cancer.xlsx")

edger_setup("prostate_cancer", counts, replicates = TRUE, filter = TRUE, gene_id = "ENSEMBL", out_path)
```

The function will organize the data into groups based on the sample labels I applied previously ("case" and "control"), filter out genes with negligible expression and calculate the expression metrics, such as the logarithm of the fold-change (logFC) and counts per million transcripts (logCPM), as well as fit a statistical generalized linear model (GLM), calculating GLM coefficients (&beta;) for each gene. The DEA then consists in perform a hypothesis test (quasi-likelihood F-test in this case), to test the null hypothesis that the coefficients are equal (or that &beta;<sub>*control*</sub> - &beta;<sub>*case*</sub> = 0). From the F-test statistics is then derived a p-value, which is adjusted by false discovery rate (FDR) to account for multiple comparisons.

After a while, the function will generate a spreadsheet with the DEA results. See below an excerpt of the spreadsheet:

![edgeR differential expression analysis in a prostate cancer dataset]({static}/images/prostate_cancer_edger_result.PNG)

Note that there are seven columns:

* `symbol` and `geneid`: transcript identifiers;
* `logfc`: the base 2 logarithm of the **fold-change** (logFC), which is how much a quantity changes between two measurements -- it is a ratio of two quantities. A logFC = 1 means that the expression of a certain gene was double in one condition than the other, a logFC = 2 means four-times higher expression, and so on. A logFC = -1 means half of the expression, a logFC = -2 means a quarter, and so on.
* `logcpm`: the logarithm of the counts per million transcripts.
* `f`: the quasi-likelihood F-test statistic.
* `pvalue` and `adjpvalue`: the quasi-likelihood F-test statistic raw and FDR-adjusted p-values, respectively.

Note that the `logfc` column is the expression in cases group relative to control group. The `edger_setup()` custom function automatically organizes data to this end.

Usually, the researcher may want to further filter these results. For example, I like to consider not only the adjusted p-value, but also check which genes presented |logFC >= 1| (note the absolute value symbols here). Thus, if a gene passes these two criteria, I usually assume that it may have biological relevance for the disease/characteristic in study.

## Conclusion

I demonstrated a custom function that uses `edgeR` package to perform differential expression analysis. Here is a summary of the requirements of the function:

* A R data frame: rows are the transcripts, columns are the samples;
* The samples must be labeled as "case" or "control" (in the column names);
* The function outputs a spreadsheet with logFC and p-values;
* The reported logFC are relative to the control group (control group is the reference);
* The result spreadsheets can be filtered as the researcher wishes.

*Subscribe to my [RSS feed](https://antoniocampos13.github.io/feeds/all.rss.xml), [Atom feed](https://antoniocampos13.github.io/feeds/all.atom.xml) or [Telegram channel](https://t.me/joinchat/AAAAAEYrNCLK80Fh1w8nAg) to keep you updated whenever I post new content.*

## References

[Differential gene expression analysis](https://www.ebi.ac.uk/training/online/course/functional-genomics-ii-common-technologies-and-data-analysis-methods/differential-gene)

[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)

[How do I explain the difference between edgeR, LIMMA, DESeq etc. to experimental Biologist/non-bioinformatician](https://www.biostars.org/p/284775/#284893)

[openxlsx package | R Documentation](https://www.rdocumentation.org/packages/openxlsx/versions/4.2.2)

[BiocManager package | R Documentation](https://www.rdocumentation.org/packages/BiocManager/versions/1.30.10)

[Bioconductor - Home](https://www.bioconductor.org/)

[AnnotationDbi](https://www.bioconductor.org/packages/release/bioc/html/AnnotationDbi.html)

[annotate package | R Documentation](https://www.rdocumentation.org/packages/annotate/versions/1.50.0)

[org.Hs.eg.db package](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html)

[EnsDb.Hsapiens.v79](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v79.html)

[Entrez Gene Records](https://www.ncbi.nlm.nih.gov/Class/MLACourse/Modules/Genes/sample_entrez_gene_record.html)

[Investigating a gene | Ensembl](https://www.ebi.ac.uk/training-beta/online/courses/ensembl-browsing-genomes/navigating-ensembl/investigating-a-gene/)

[HUGO Gene Nomenclature Committee](https://www.genenames.org/)

[Human housekeeping genes, revisited](https://www.cell.com/trends/genetics/fulltext/S0168-9525(13)00089-9?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0168952513000899%3Fshowall%3Dtrue)

[Human housekeeping genes, revisited - Supplementary material](https://www.tau.ac.il/~elieis/HKG/)
