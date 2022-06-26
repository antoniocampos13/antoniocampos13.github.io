---
Title: Opening files of size larger than RAM with pandas
Status: published
Date: 2022-06-27 10:00
Author: Antonio Victor Campos Coelho
Categories: Python
Tags: pandas, Genomics, Bioinformatics
---

## Introduction

Dealing with big files is a routine for everyone working in genomics. FASTQ, VCF, BAM, and GTF/GFF3 files, to name a few, can range from some hundreds of megabytes to several gigabytes in size. Usually, we can use cloud services to configure computing instances with a lot of RAM, but we may use some ways to read and manipulate large-than-RAM files in our personal/work machines.

This post will demonstrate how to work with big tabular data using the `chunksize` option with `pandas`. You can find the code in my [portfolio](https://github.com/antoniocampos13/portfolio/tree/master/Python/2022_06_27_Opening_files_of_size_larger_than_RAM_with_pandas).

### Chunking: divide and conquer

"Chunking" means splitting the big file into chunks (partitions) so the Python session can work with each part separately, meaning it would not need to hold the big data in memory all at once. Keep in mind that not every problem can be solved by chunking. Therefore, if your goal does NOT involve coordination between chunks, such as some filtering and little edition, chunking could help. However, if your task is more complicated than this, other modules such as [`dask`](https://www.dask.org/) are the better option. The panda's documentation has an excellent [chapter](https://`pandas`.pydata.org/docs/user_guide/scale.html#) explaining ways to go when scaling to large datasets.

## The input: the human reference genome GTF/GFF3 file

I downloaded the human reference genome GTF/GFF3 file `GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz` file at the [NCBI FTP server](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) containing files preformatted for use in Bioinformatic analysis pipelines. Next, I extracted the contents of the file:

```bash
gunzip GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
```

The extracted file has a size of about 1 GB. Decently big for my demonstrational purposes.

But what is a GTF/GFF3 file? The Gene Transfer Format or General Feature Format is a tab-delimited text file format. Bioinformaticians use it to describe genomic features such as genes, exons, introns, putative protein-coding sequences (CDS), transcription factor binding sites, etc. The first two versions (GTF and GFF2) had deficiencies, and GFF3 was developed to address them. You can read more about GFF3 on [Lincoln Stein's GitHub page](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md).

Every GTF/GFF3 file has nine fields (columns). A dot `.` represents missing or null data. The nine columns are:

1. `seqid`: the name of the sequence where the feature is located. For example, a chromosome or contig;
2. `source`: the program or organization, laboratory, etc. that generated the information regarding the feature;
3. `type`: qualifiers like "gene", "exon", "CDS". Features can have children: for example, the exons of a gene refer to its gene (their parent). Ideally, all the children features must follow their parents after the parents' initial definition in the file.
4. `start`: the base position in the sequence where the feature starts. It has a 1-base offset, in contrast to the BED format, which is 0-offset.
5. `end`: the base position in the sequence where the feature ends.
6. `score`: numeric value representing the quality of the sequence.
7. `strand`: indicates the strand of the feature: `+` (the sense strand is the default 5'-3' representation of the feature), `-` (the sense strand is the reverse complement strand of the sequence representation), or `.` (undetermined).
8. `phase`: used to indicate the reading frame of the features that are CDS. Can be `0`, `1`, `2` or `.` (undetermined).
9. `attributes`: all other information relevant for describind the feature.

I will open this file using `pandas` and keep just the exons of all annotated human genes.

## Chunking with `pandas`

I create a python script and import the `pandas` module:

```python
import `pandas` as pd
```

Next, I define some variables to store the path of the GTF/GFF3 file, the name of the output table, and the chunk size:

```python
TABLE = "GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff"
EDITED_TABLE = "human_gene_exons_GRCh38.gff"
CHUNKSIZE=20000000
```

The chunk size must be an integer because it represents the number of lines each chunk will have. In the example above, I will tell `pandas` to partition the file into parts containing two million lines until it finishes processing the whole dataset.

Next, I define a list with the column names:

```python
column_names = [
    "seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes"
]
```

Since the file is tab-delimited, I can use the `pd.read_table()` function, passing the `column_names` list as the value for the `names` argument and the `chunksize` as well:

```python
chunks = pd.read_table(TABLE, names=column_names, chunksize=CHUNKSIZE)

# Or:
chunks = pd.read_csv(TABLE, names=column_names, chunksize=CHUNKSIZE, sep="\t")
```

If I had not decompressed the file, I could also use the argument `compression="infer"` so `pandas` would decompress it on-the-fly.

The `chunksize` argument makes the `pd.read_table()` return an **iterator** object (`TextFileReader`). What is an iterator? In Python, an iterator is an object we can traverse through all its values. Python lists, dictionaries, and tuples are all Pythonic iterators. In our case, each element of this iterator is a `pd.Dataframe` instead.

Therefore, the `for` loop I wrote will perform the same action on every `pd.DataFrame` in the iterator.

## Filtering each chunk

These are the steps I will perform with each chunk:

1. Filter for rows with "exon" values in the `type` column;
2. Drop (remove) all columns except for `seqid`, `type`, `start`, and `end`;
3. Save the edited chunk directly to disk by appending the chunk to a tab-delimited file.

Below is the loop code:

```python
for chunk in chunks:
    # Step 1
    temp_df = chunk[chunk["type"] == "exon"]

    # Step 2
    temp_df = temp_df.drop(["source", "score", "strand", "phase", "attributes"])

    # Step 3
    temp_df.to_csv(EDITED_TABLE, header=False, mode="a", sep="\t", index=False)
```

I should explain step 3 in more detail. Since the `for` loop will dump each part immediately to disk after finishing my edits, Python's RAM usage will be more or less constant during the file processing. You may see other tutorials appending the data frames to a list and concatenating them into a final `pd.DataFrame`, but in my opinion, this kind of defeats the purpose of chunking since Python will have to hold everything in memory, risking RAM overuse and killing the process.

Let me explain the `pd.to_csv()` arguments:

- `EDITED_TABLE`: the file output name/path;
- `header=False`: do not output the column names to file;
- `mode="a"`: append each chunk on the output file;
- `sep="\t"`: write tab-delimited columns on the output file;
- `index=False`: do not output index column. Since I did not set the index, it would print the row numbers, which would be undesirable (it would violate GTF/GFF3 format specifications).

Observe that I defined the column names to make filtering/editing easier. The GTF/GFF3 format specifications do not require the header names to be present in the file. Therefore, I removed them during Step 3.

## Conclusion

In this post, I demonstrated one way of dealing with big files by chunking with `pandas`.

## Appendix

I already have written about `dask` on my machine learning tutorials. See part 1 [here](https://antoniocampos13.github.io/machine-learning-with-python-supervised-classification-of-tcga-prostate-cancer-data-part-1-making-features-datasets.html) and part 2 [here](https://antoniocampos13.github.io/machine-learning-with-python-supervised-classification-of-tcga-prostate-cancer-data-part-2-making-a-model.html). Check [`pandas`' documentation](https://pandas.pydata.org/docs/index.html) and [dask's](https://docs.dask.org/en/stable/) as well for more information.

*Subscribe to my [RSS feed](https://antoniocampos13.github.io/feeds/all.rss.xml), [Atom feed](https://antoniocampos13.github.io/feeds/all.atom.xml) or [Telegram channel](https://t.me/joinchat/AAAAAEYrNCLK80Fh1w8nAg) to keep you updated whenever I post new content.*

## References

[Dask | Scale the Python tools you love](https://www.dask.org/)

[Scaling to large datasets &#8212; pandas 1.4.3 documentation](https://pandas.pydata.org/docs/user_guide/scale.html#)

[Index of /genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/)

[Specifications/gff3.md at master · The-Sequence-Ontology/Specifications](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)

[Machine Learning with Python: Supervised Classification of TCGA Prostate Cancer Data (Part 1 - Making Features Datasets)](https://antoniocampos13.github.io/machine-learning-with-python-supervised-classification-of-tcga-prostate-cancer-data-part-1-making-features-datasets.html)

[Machine Learning with Python: Supervised Classification of TCGA Prostate Cancer Data (Part 2 - Making a Model)](https://antoniocampos13.github.io/machine-learning-with-python-supervised-classification-of-tcga-prostate-cancer-data-part-2-making-a-model.html)
