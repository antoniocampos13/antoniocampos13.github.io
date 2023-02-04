---
Title: Parsing the ClinVar XML file with pandas
Status: published
Date: 2023-02-04 16:25
Author: Antonio Victor Campos Coelho
Categories: Python
Tags: pandas, ClinVar, genomics, variants
---

## Introduction

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/intro/) is one of the USA's National Center for Biotechnology Information (NCBI) databases. ClinVar archives reports of relationships among human genetic variants and phenotypes (usually genetic disorders). Any organization, such as a laboratory, hospital, clinic etc can submit data to ClinVar. The core idea of ClinVar is aggregate evidence for the clinical significance of any genetic variant concerning any disorder. Over 2,400 organizations contributed more than 2 million 600 thousand [records to ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/submitters/), representing more than 1 million 600 thousand unique variants.

Anyone can freely search ClinVar through their [website](https://www.ncbi.nlm.nih.gov/clinvar/), using gene symbols, genomic coordinates, [HGVS expressions](https://varnomen.hgvs.org/), phenotypes, and more. If you want to perform a few queries, the online search tool does a good job. However, if you are pursuing more complex scientific questions, or are intending to download batches of data, the search tool will not suffice. Other NCBI databases can be queried via the command line with the [Entrez Direct (EDirect) utilities](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (in a [previous post](https://antoniocampos13.github.io/fastq-to-annotation-part-2.html) I mention how to work with the EDirect utilities). Unfortunately, ClinVar does not currently support a batch query interface via EDirect utilities.

However, ClinVar provides [other approaches](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/) for the access and use of their data. One of these approaches is the provisioning of the complete public data set in the form of an [XML file stored at the ClinVar FTP server](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/). The `ClinVarFullRelease` XML file is updated weekly, and every release happening on the first Thursday of the month is archived.

## Parsing the ClinVar XML file

Recently, I started assisting my team in uploading variant/phenotype interpretations to ClinVar. I wanted to find a way to gather all our submissions into a spreadsheet so every team member could easily check whenever necessary. Thus, I downloaded the full ClinVar release XML file and tried to parse it with the XML-handling [`ElementTree` module](https://docs.python.org/3/library/xml.etree.elementtree.html). However, I had limited success. I could extract some information, but the output did not turn out exactly the way I was intending, so I set out to find working alternatives.

Eventually, I found out that the [`pandas` module](https://pandas.pydata.org/) has a method to convert XML-stored data into traditional data frames. Moreover, since September 2022, their `read_xml()` function supports large XML files via the `iterparse` argument (read an excerpt of the release note [here](https://pandas.pydata.org/docs/whatsnew/v1.5.0.html#read-xml-now-supports-large-xml-using-iterparse)).

The function documentation states that the `iterparse` argument is a memory-efficient method for handling big XML files without storing all data elements within memory. This was exactly my case, so I tried the `read_xml()` function &mdash; it worked quite well!

I wrote a small script that you can use to parse the ClinVar XML file. Of course, when you get acquainted with the `read_xml()`, you may use it for parsing any other XML you wish. I used an [AWS EC2 instance](https://aws.amazon.com/ec2/?nc1=h_ls) with 90 GB RAM while working on this tutorial. I did not try to process the ClinVar file in less powerful systems. *Try at your own risk*.

I uploaded the script (named `clinvar_pandas_xml_parser.py`) to [the corresponding folder on my portfolio](https://github.com/antoniocampos13/portfolio/tree/master/Python/2023_02_04_Parsing_ClinVar_XML_with_pandas).

### Downloading the full ClinVar release XML file

You can download the latest `.gz`-compressed XML files via the following links (**WARNING:** the release files are HUGE):

[Weekly release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarFullRelease_00-latest_weekly.xml.gz)

[Monthly release](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz)

Go to a convenient directory on your system and download one of the files above. I downloaded the most recent monthly release and decompressed it soon after:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz

gunzip ClinVarFullRelease_00-latest.xml.gz
```

If you want to check the file integrity, compare your checksum against the corresponding ClinVar-provided checksum file:

[Weekly release (MD5 checksum file)](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/weekly_release/ClinVarFullRelease_00-latest_weekly.xml.gz.md5)

[Monthly release (MD5 checksum file)](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarFullRelease_00-latest.xml.gz.md5)

### Installing modules

The `iterparse` argument in the `read_xml()` function was introduced in `pandas` version 1.5.0 and requires the `lxml` or `ElementTree` modules to work. In this tutorial, I will use `lxml` (the default). Therefore, install the necessary modules via `pip` or `conda` (see my [previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html) on how to configure `conda` virtual environments in a Unix system). For example:

```bash
conda activate env_name
conda install -c conda-forge pandas=1.5.0 lxml
```

### Running the `clinvar_pandas_xml_parser.py` script

Finally, let's walk through the script. First, I import the `pandas` module:

```python
import pandas as pd
```

Next, I saved the XML file path into the `xml_file_path` object:

```python
xml_file_path = "ClinVarFullRelease_00-latest.xml"
```

Then, I investigated the XML using `grep` commands to match specific strings of interest to get a feel of how the XML file was structured. I am sure that are better ways to assess the XML elements structure, but I am not an expert in XML files.

Through my investigation of the file, I concluded that the `ClinVarAssertion` elements within the XML structure contained all information I was needing at the moment. Thus, I created a Python dictionary object named `iterparse_dict` with the string "`ClinVarAssertion`" as a *key* and a Python list as its corresponding *value*:

```python
iterparse_dict = {"ClinVarAssertion": []}
```

The *key* represents the parent XML node tag. The *value* is a list containing all child or grandchild nodes, tags, or attributes at any node level inside the main XML node &mdash; simple as that. I chose the following:

```python
iterparse_dict = {
    "ClinVarAssertion": [
        "ID",
        "SubmissionName",
        "localKey",
        "submittedAssembly",
        "submitter",
        "submitterDate",
        "Acc",
        "RecordStatus",
        "OrgID",
        "DateCreated",
        "DateUpdated",
        "Version",
        "DateLastEvaluated",
        "Description",
        "ReviewStatus",
        "Comment"
    ]
}
```

Then, I passed the `iterparse_dict` as the value for the `iterparse` argument of the `read_xml()` function and stored the output as the `df` object &mdash; a `pandas.DataFrame`. The columns of the data frame will correspond to the information stored at each `ClinVarAssertion` tag, attributes

```python
df = pd.read_xml(xml_file_path, iterparse=iterparse_dict)
```

After some time, the function returns a data frame that you can further filter to search for information. For now, I saved the data frame as a pickled object:

```python
df.to_pickle("pandas_parsed.pkl")
```

At any moment, I can restore the data frame through `pandas` as well (**REMEMBER:** Loading pickled data received from untrusted sources can be unsafe):

```python
df = pd.read_pickle("pandas_parsed.pkl")
```

## Conclusion

In this post, I demonstrated one way of exploring the full release of the ClinVar database, through an up-to-date `pandas` method that can deal with big XML files.

## Appendix

A brief explanation of what each *value* in the `iterparse_dict` dictionary object represents:

- `ID`: A unique numeric id representing a submission (a single submission usually contains many variant/phenotype interpretations).
- `SubmissionName`: A unique string representing a submission.
- `localKey`: The HGVS expression representing each variant within a single submission.
- `submittedAssembly`: The assembly (genome reference build) that was used for variant calling, annotation and localization. Usually is "GRCh37" or "GRCh38".
- `submitter`: The organization that was responsible for the submission.
- `submitterDate`: The date when the submission was uploaded to ClinVar.
- `Acc`: A ClinVar identifier string. As stated in the [ClinVar identifiers documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/): "Accession numbers in ClinVar have the pattern of 3 letters and 9 numerals. The letters are either SCV (think of it as Submitted record in ClinVar), RCV (Reference ClinVar record) or VCV (Variation ClinVar record)."
- `RecordStatus`: The status of the record, whether current, deleted or secondary (merged).
- `OrgID`: A unique numeric identifier for each organization that was responsible for the submission.
- `DateCreated`: The date when the submission was accepted and integrated into the database.
- `DateUpdated`: The date when the submitter updated the record.
- `Version`: The version assigned to a record. As stated in the [ClinVar identifiers documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/): "The version number is incremented when a submitter updates a record or when the contents of a reference or variation record change because of addition to, updates of, or deletion of the SCV accessions on which it is based."
- `DateLastEvaluated`: The date when the organization evaluated the clinical significance of any given variant in the context of any given phenotype.
- `Description`: A description of the clinical significance of any given variant in the context of any given phenotype, such as "Pathogenic", "Likely pathogenic", "Benign", etc.
- `ReviewStatus`: As stated in the [ClinVar review status documentation](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/): "The level of review supporting the assertion of clinical significance for the variation."
- `Comment`: Any additional (free-text) comments the organization that was responsible for the submission provided regarding the interpretation of any given variant in the context of any given phenotype.

*Subscribe to my [RSS feed](https://antoniocampos13.github.io/feeds/all.rss.xml), [Atom feed](https://antoniocampos13.github.io/feeds/all.atom.xml) or [Telegram channel](https://t.me/joinchat/AAAAAEYrNCLK80Fh1w8nAg) to keep you updated whenever I post new content.*

## References

[ClinVar | Documentation | What is ClinVar?](https://www.ncbi.nlm.nih.gov/clinvar/intro/)

[ClinVar | Documentation | Submitters](https://www.ncbi.nlm.nih.gov/clinvar/submitters/)

[ClinVar | Search Tool](https://www.ncbi.nlm.nih.gov/clinvar/)

[Sequence Variant Nomenclature](https://varnomen.hgvs.org/)

[Entrez Direct: E-utilities on the Unix Command Line - Entrez Programming Utilities Help - NCBI Bookshelf](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

[FASTQ to Annotation (Part 2)](https://antoniocampos13.github.io/fastq-to-annotation-part-2.html)

[Accessing and using data in ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)

[ClinVar | FTP server | Index of /pub/clinvar/xml](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/)

[xml.etree.ElementTree — The ElementTree XML API](https://docs.python.org/3/library/xml.etree.elementtree.html)

[pandas - Python Data Analysis Library](https://pandas.pydata.org/)

[pandas | Documentation | What’s new in 1.5.0 (September 19, 2022)](https://pandas.pydata.org/docs/whatsnew/v1.5.0.html#read-xml-now-supports-large-xml-using-iterparse)

[Secure and resizable cloud compute – Amazon EC2 – Amazon Web Services](https://aws.amazon.com/ec2/?nc1=h_ls)

[Setting Up Your Unix Computer for Bioinformatics Analysis](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html)

[ClinVar | Documentation | Identifiers in ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/docs/identifiers/)

[ClinVar | Documentation | Review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
