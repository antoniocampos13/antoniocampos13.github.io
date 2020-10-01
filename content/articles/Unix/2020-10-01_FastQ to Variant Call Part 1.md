---
Title: FASTQ to Variant Call (Part 1)
Date: 2020-10-01 18:00
Author: Antonio Victor Campos Coelho
Tags: Bioinformatics, genomic variation, entrez-direct, EDirect
---

## Introduction

On my [previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html#setting-up-your-unix-computer-for-bioinformatics-analysis), I showed how to configure an Ubuntu system to install Bioinformatics programs.

Now, using the environment I created, I will demonstrate a bash script that takes next generation sequencing (NGS) raw reads from human whole genome sequencing as input and produces variant annotation as output. Variant annotation is the process of identifying genetic variants in some genomic DNA sample, and assess, for example, if any of the found variants have any effect on phenotype, such as increased susceptibility to certain diseases.

This demonstration will be separated in parts. In the first part, I will show how to search for NGS projects deposited in [National Center for Biotechnology Information (NCBI) databases](https://www.ncbi.nlm.nih.gov/) from which we can download sequencing reads later.

## Using NCBI's entrez-direct (EDirect) to retrieve FASTQ files

I open my Unix terminal and activate the `bioenv` environment:

```bash
conda activate bioenv
```

Now I use the `EDirect` `esearch` command to search NCBI's databases. We must provide a database using the flag `-db`. Let's first check the available databases. If I type:

```bash
esearch -db
```

Here is an excerpt of the output showing the available databases:

```bash
# COLOCAR OUTPUT
```

I will search the `biproject` database because it contains metadata from projects dealing with high-throughput genome sequencing, transcriptome expression analysis, and several other designs. We must use the `-query` flag to provide keywords for search. In this example, I will search for studies dealing with **vorinostat**, a medicine that is have been used in experimental HIV-1 latency reversal, or "shock-and-kill" treatments.

Remember to use single quotes ('') enclosing the query, especially if it has several words.

```bash
# It is just the beginning... (1/4)
esearch -db bioproject -query 'vorinostat'
```

The output is only a summary of the number of results retrieved. Thus, I will add more commands to retrieve the actual query results. I will pipe, i.e. transfer, the results of the query to the another command -- `efetch` -- that will do this work for me. This is the pipe symbol: (`|`)

```bash
# ... not there yet ... (2/4)
esearch -db bioproject -query 'vorinostat' | efetch -format native -mode xml
```

 The output is in `XML` format, and it is unfortunately not very much human-readable. Thus, I will once again pipe the results, this time to `xtract` command. As its name implies, it extracts information from the `XML` and formats into a tab-separated format that is easier to understand. We must input the flag `-pattern` with the part of the `XML` files that contains the desired information, which are `elements`. In this example, I will search inside the `DocumentSummary` for `ArchiveID@accession` (project unique accession number), `ID` (an auxiliary ID code to search for samples of said project), `Title`(the title of the project),  `Description` (normally an abstract of the project) and `Reference` (a list of project-related papers in PubMed ids -- PMIDs, if available). Note that we are separating each argument with spaces, no quotes are necessary in this command.

 ```bash
# ... almost there ... (3/4)
esearch -db bioproject -query 'vorinostat' | efetch -format native -mode xml | xtract -pattern DocumentSummary -element ArchiveID@accession ID Title Description Reference
 ```

The output is displayed in the terminal. Lastly, I will add a final command to transfer to a local text file `vorinostat_projects.txt` that will be saved in the current working directory.  Note that if you have a identically-named file in the working directory, it will be overwritten, so be careful.

```bash
# Finally there! (4/4)
esearch -db bioproject -query 'vorinostat' | efetch -format native -mode xml | xtract -pattern DocumentSummary -element ArchiveID@accession ID Reference Title Description > vorinostat_projects.txt
```

## Refining the search

Thus, we have a basic command to search NCBI databases via `EDirect`. We can create more elaborate queries by adding other keywords and filtering results. NCBI's search engines have several parameters. I advise you go to any advanced search page on the NCBI website to look for the available parameters.

Using `bioproject` as example again:

Click on *Advanced* to go the query constructor