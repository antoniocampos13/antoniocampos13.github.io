---
Title: Genomic Analysis With Hail
Status: published
Date: 2021-07-09 16:45
Author: Antonio Victor Campos Coelho
Categories: Python
Tags: Bioinformatics, Genomics, Hail
---

## Introduction

Hello, long time no see! Since I lasted posted, many things happened. Since March I have been working as Post-Doc Researcher, hired by the [Hospital Israelita Albert Einstein (HIAE, São Paulo, Brazil)](https://www.einstein.br/Pages/Home.aspx) to work for the Projeto Genomas Raros ("Rare Genomes Project", GRAR from here on), a public-private partnership between HIAE and the Brazilian Health Ministry to further the implementation of genomic analysis into the Brazilian public healthcare system (SUS), with the intention of improve diagnostic rates of rare diseases in Brazil. Since 2020, thousands of genomes of Brazilian patients with suspected rare diseases have been sequenced, and many more will come in the next two years.

Thus, I have been tasked to develop/adapt analysis pipelines to handle whole-genome data at large scales compatible with the scope of GRAR. My team asked me to explore [Hail, a Python/Spark framework for scalable genomic analysis](https://hail.is/). It has been developed at [Broad Institute](https://www.broadinstitute.org/) and was used to generate the [Genome Aggregation Database (gnomAD)](https://gnomAD.broadinstitute.org/).

In this post I will share some use cases of this tool full of potential. Please notice that it is not intended to substitute the official documentation of Hail, GATK and other software used here. Just consider it as a demonstration of Hail use cases with commentaries. Also, notice two important things: first, Hail implements "lazy evaluation"; new users may think a specific command run blazingly fast, but in reality, Hail just mapped the execution order of functions needed and only will compute anything when necessary, such as saving results to disk and printing the first few rows of a dataset to Python's standard output stream. Second, Hail need a lot of RAM and CPU to work properly with big datasets. Hail is intended to be used in cloud/cluster computing environments, but a computer with relatively good hardware configuration can run small datasets.

## Installing software

Prepare your Unix computing environment by installing [Hail](https://hail.is/#install), [gnomAD utilities for Hail](https://pypi.org/project/gnomAD/) and [Genome Analysis Toolkit version 4 (GATK4)](https://gatk.broadinstitute.org/hc/en-us). Tools for manipulating VCF files are essential too, such as [bcftools](http://samtools.github.io/bcftools/bcftools.html). I advise installing the tools into a virtual environment for convenience (see [my previous post](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html) for some pointers on how to use conda environments).

## Preparing the multi-sample VCF input

Since the GRAR data is coming from whole-genome sequencing (WGS), we have been generating one [genomic variant call format (gVCF) file](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format) per participant through [Illumina's DRAGEN sequencing read analysis platform](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html). A gVCF file has all the characteristics of a [VCF file](https://samtools.github.io/hts-specs/VCFv4.2.pdf), the difference being that the gVCF files have information for all sites in the genome, which give more precision during the integration of variant calls coming from several samples.

By integration I mean **combining** several gVCF files into a single multi-sample VCF file so we can perform analysis (calculate allelic frequencies, assess genotyping quality, and so on) with the **whole cohort**. We are currently using GATK's GenomicsDBImport tool. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport) for a tutorial of how to use it. Briefly, we intend to create a GenomicsDB object and we will update it regularly by appending new gVCFs as they become available, until the last participant is recruited and we have their genome sequenced. When this moment comes, we will extract a multi-sample VCF with GATK's `GenotypeGVCFs` tool:

```bash
gatk GenotypeGVCFs -R $REF -V $DBPATH -G StandardAnnotation -O cohort.vcf
```

Where `$REF` and `$DBPATH` are the paths of the genome reference (the same used during the variant call process) and the GenomicsDB, respectively (which should have a `gendb://` prefix as noted in the GenomicsDBImport tutorial, something like `gendb://my_database`). The `-O` flag indicates the output name. If you wish to compress the VCF file, you may use `bcftools` or `bgzip`. Remember to index the compressed file:

```bash
bcftools view cohort.vcf -Oz -o cohort.vcf.gz
bcftools index cohort.vcf.gz

# or
bgzip -@ 4 cohort.vcf
tabix cohort.vcf.gz
```

The next step is to perform GATK's Variant Quality Score Recalibration (VQSR) in the output VCF with GATK's `VariantRecalibrator`.  See their original tutorial [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering). To put it simply, VQSR works by comparing the detected variants with high-confidence variant sites observed by several consortia (HapMap, 1000 Genomes etc.) and applies a filter deeming the variant a true positive (i.e. the observed variation is a true biological event) or a false positive (i.e. the observed variation is not real, it is in a fact sequencing artifact). To this end, I downloaded high-confidence datasets to perform the VQSR. The datasets can be found at the [GATK Google Cloud Storage (Resource Bundle)](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false) and can be downloaded with [`gsutil` application](https://cloud.google.com/storage/docs/gsutil).

I will post a script with slight modification in the steps of the VQSR tutorial in my [portfolio](https://github.com/antoniocampos13/portfolio/tree/master/Python/2021_07_07_Genomic_Analysis_With_Hail) (I was having errors so I noticed that I had to put spaces after the `-resource` flags of the `VariantRecalibrator` command). Briefly, the steps are:

1. Filtering samples with excess of heterozygotes (recommended when working with thousands of samples);
2. Make a sites-only VCF;
3. Recalibration step (separately by indels/mixed and SNP loci);
4. Apply the recalibration filters.

The output of the VQSR is the `snp.recalibrated.vcf.gz` file (despite the name, the indels/mixed variants in the file have been recalibrated as well). We can now import the dataset into Hail.

## Initiating Hail

Hail's frontend is written in Python. Simply importing the Hail module is not sufficient. We must initiate it so it can communicate with Spark. I created the `hail_demo_init.py` file as an "init template" that can be reused between Hail scripts/sessions:

```python
# hail_demo_init.py
import hail as hl

DEFAULT_REF = "GRCh38"

hl.init(
    idempotent=True,
    quiet=True,
    skip_logging_configuration=True,
    default_reference=DEFAULT_REF,
    spark_conf={
        "spark.executor.cores": "4",
        "spark.driver.memory": "16g",
    },
)
```

Notice the `DEFAULT_REF` variable: it establishes that I intend to use the genomic coordinates considering the Genome Reference Consortium Human Reference 38 (GRCh38), the most recent genome version being accepted by the human genome community. Otherwise Hail would use GRCh37 as default. The `spark_conf` argument where I setup Spark so it uses four CPU cores and 16 GB of RAM. Change these values as appropriate.

## Importing the VCF input

I am now ready to import `snp.recalibrated.vcf.gz` into the Hail session and convert it to Hail's `MatrixTable` object and write it to the disk. In simple terms, a `MatrixTable` is a representation of a VCF file that is amenable to be manipulated by Spark. Read the [Hail Docs](https://hail.is/docs/0.2/index.html) for more details. I initiate Hail and then import and convert the dataset with two chained steps (`hl.import_vcf().write()`):

```python
# hail_demo_import_vcf.py
import hail as hl
from hail_demo_init import DEFAULT_REF

hl.import_vcf("snp.recalibrated.vcf.gz",
force_bgz=True,
reference_genome=DEFAULT_REF,
array_elements_required=False).write("recalibrated.mt",overwrite=True)
```

I can now proceed to prepare the dataset for sample and variant quality control (QC).

## Preparing for Sample QC

Before reading the created matrix table I create some "magic numbers" variables to hold some quality metrics thresholds for sample QC later:

```python
# hail_demo_sample_qc.py
import hail as hl
from hail_demo_init import DEFAULT_REF

from gnomAD.utils.filtering import filter_to_autosomes
from gnomAD.utils.annotations import add_variant_type

# Magic numbers
CALL_RATE = 0.90  # Hail team value = 0.97. gnomAD value = 0.99
RELATEDNESS = 0.088  # Hail team value. gnomAD value = 0.08838835
READ_DEPTH = 20  # Hail team value.
FREEMIX_CONTAMINATION = 0.05  # gnomAD value.
CHIMERIC_READS = 0.05  # gnomAD value.
MEDIAN_LENGTH = 250  # gnomAD value.
```

The comment in each line contains the value used by Hail team or gnomAD team. Since this is a demonstration, I may have used different values. Change these values as you feel appropriate as well.

I will now read the matrix table from disk and assign it to the `mt` object:

```python
mt = hl.read_matrix_table("recalibrated.mt")
```

To check the first few sample ids, I use the `show()` method (I explain the `s` here later):

```python
mt.s.show()
```

To check the matrix table structure, I use the `describe()` method:

```python
mt.describe()
```

The output is as follows:

```text
----------------------------------------
Global fields:
    None
----------------------------------------
Column fields:
    's': str
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'alleles': array<str>
    'rsid': str
    'qual': float64
    'filters': set<str>
    'info': struct {
        AC: array<int32>, 
        AF: array<float64>, 
        AN: int32, 
        BaseQRankSum: float64, 
        DB: bool, 
        DP: int32, 
        END: int32, 
        ExcessHet: float64, 
        FS: float64, 
        FractionInformativeReads: float64, 
        InbreedingCoeff: float64, 
        LOD: float64, 
        MLEAC: array<int32>, 
        MLEAF: array<float64>, 
        MQ: float64, 
        MQRankSum: float64, 
        NEGATIVE_TRAIN_SITE: bool, 
        POSITIVE_TRAIN_SITE: bool, 
        QD: float64, 
        R2_5P_bias: float64, 
        ReadPosRankSum: float64, 
        SOR: float64, 
        VQSLOD: float64, 
        culprit: str
    }
----------------------------------------
Entry fields:
    'AD': array<int32>
    'AF': array<float64>
    'DP': int32
    'F1R2': array<int32>
    'F2R1': array<int32>
    'GP': array<float64>
    'GQ': int32
    'GT': call
    'ICNT': array<int32>
    'MB': array<int32>
    'MIN_DP': int32
    'PL': array<int32>
    'PRI': array<float64>
    'PS': int32
    'RGQ': int32
    'SB': array<int32>
    'SPL': array<int32>
    'SQ': float64
----------------------------------------
Column key: ['s']
Row key: ['locus', 'alleles']
----------------------------------------
```

We can see that a Hail `MatrixTable` objects has four types of information "compartments":

- Global fields: information values that are identical for every row
- Column fields: sample-level information (id, phenotype, sex etc.)
- Row fields: variant-level information (locus, alleles, type of variant, allele frequency etc.)
- Entry fields: variant-by-sample-level (genotype, genotype quality, etc.)

Also notice that the `MatrixTable` is keyed by column (`s`: **s**ample id field) and rows (locus and allele row fields) allowing us to perform SQL-style table joins. For each field type you can see the name of each field. The info field contains the INFO field from the input VCF. The entries fields contain the values listed into the FORMAT VCF field. For more details check the [`MatrixTable` Hail Docs](https://hail.is/docs/0.2/tutorials/07-matrixtable.html).

Now I will modify an entry field and create other row fields that I will need to use later during variant QC. The commands to create/modify fields is `hl.annotate_cols()`, `hl.annotate_rows()` or `hl.annotate_entries()` depending on the field type (columns, rows, entries, respectively). Notice that in one of them I used a gnomAD function.

```python
# Mixture of non-empty with empty PL fields causes problems with sample QC for some reason; setting field to all empty
mt = mt.annotate_entries(PL=hl.missing(mt.PL.dtype))

# Add variant-level annotations necessary for variant QC later
## Annotate variants in one of the categories: SNV, multi-SNV, indel, multi-indel, mixed
mt = mt.annotate_rows(**add_variant_type(mt.alleles)) # gnomAD function

## Number of alleles at the site
mt = mt.annotate_rows(n_alleles = hl.len(mt.alleles))

## Mixed sites (SNVs and indels present at the site)
mt = mt.annotate_rows(mixed_site = hl.if_else(mt.variant_type == "mixed", True, False))

## Spanning deletions
mt = mt.annotate_rows(spanning_deletion=hl.any(lambda a: a == "*", mt.alleles))
```

To check the dimensions of the `MatrixTable`, I can use the following commands:

```python
# Number of Rows, Columns
mt.count()

# Number of Columns
mt.count_cols()
```

To see a variant breakdown (types, number of variants per chromosome and several other information), there is the following command:

```python
hl.summarize_variants(mt)
```

A recommended step for further downstream analyses is to split multiallelic variants into biallelic configuration:

```python
mt = hl.split_multi_hts(mt)
```

To remove any monomorphic (invariant) loci I use:

```python
mt = mt.filter_rows(mt.n_alleles > 1)
```

Try the `hl.summarize_variants()` to check how the numbers changed after splitting.

Up until this moment, I have only genetics-related information into the `MatrixTable`. I can instruct Hail to import TSV/CSV file with sample-level (column) annotations with the `hl.import_table()` command and then associate each observation by keying the `s` field with `hl.annotate_cols()`.

```python
sa = hl.import_table("sample_info.txt", impute=True, key="s")

mt = mt.annotate_cols(sample_info=sa[mt.s])
```

The `sample_info.txt` example file has the following format:

```text
s sex phenotype mean_coverage chimeric_reads contamination median_length
sample1 XX case 25 0.02 0.00 300
sample2 XY case 22 0.01 0.00 250
sample3 XY control 30 0.03 0.00 265
sample4 XX control 29 0.01 0.00 250
```

In this mock sample information file I included the sex karyotype, a fictitious phenotype and sequencing metrics to support the sample QC procedure. As you may have guessed, any number and type (numeric, string, float, integer, Boolean) of important sample metadata columns can be included.

Whenever we annotate columns, rows or entries in Hail, we must provide the name of the new field. In this case is `sample_info`. So if I wanted to manipulate, say the sex karyotype field, I would refer the field name this way:

```python
mt.sample_info.sex
```

This is because the information in the file will be assigned to a dictionary-like field named `sample_info` within the `MatrixTable`; the `sex`, `phenotype`, etc. fields are *inside* `sample_info`.

## Performing Sample QC

### Plotting quality metrics

Hail has a very convenient function to calculate descriptive statistics of fields in a `MatrixTable`:

```python
mt = hl.sample_qc(mt)
```

This function will calculate overall call rate per sample, mean read depth (coverage) and other [things](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc). We can plot the calculated fields to assess call rate and coverage (read depth) across samples. The code below will output a HTML plot to disk (make sure you have the [bokeh](https://docs.bokeh.org/en/latest/index.html) Python module installed).

```python
from bokeh.embed import file_html
from bokeh.resources import CDN

p = hl.plot.scatter(
    x=mt.sample_qc.dp_stats.mean,
    y=mt.sample_qc.call_rate,
    xlabel="Mean DP",
    ylabel="Call Rate",
    hover_fields={"ID": mt.s},
    size=8,
)

html = file_html(p, CDN, "Chart")

with open("Call Rate by Mean DP.html", "w") as f:
    f.write(html)
```

### Filter samples by quality thresholds

I can finally filter out samples that do not meet the quality thresholds I established before. For example, the lines below will keep only samples that meet overall call rate and mean coverage quality criteria:

```python
mt = mt.filter_cols(mt.sample_qc.call_rate >= CALL_RATE)

mt = mt.filter_cols(mt.sample_qc.dp_stats.mean >= READ_DEPTH)
```

If you have other quality criteria, you could filter in a similar way by correctly referring the field name (remember that dictionary-like fields contain other fields, as is the case of `sample_info` I mentioned earlier and `sample_qc` above) and a logical expression (which must follow the Python syntax for equalities and inequalities).

### Principal component analysis (PCA) to filter related samples

To ensure that each sample in the cohort is unrelated to any other sample, we can run a principal component analysis (PCA) to evidence samples with kinship too higher according to our chosen threshold. The PCA will only work with autosome (diploid) biallelic variants. The code below will filter our `MatrixTable` keeping only variants meeting these conditions:

```python
for_pca = filter_to_autosomes(mt) # gnomAD function. Will keep variants in chromosomes 1 to 22 only.
for_pca = for_pca.filter_rows(for_pca.n_alleles == 2) # will remove sites that were multi-allelic before splitting as well to ensure "pure" biallelic sites
```

Then, I determine the sample number to calculate `k`, the number of principal components.

```python
sample_num = for_pca.cols().count()
```

Next, using the [`hl.hwe_normalized_pca()` function](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca) I calculate the principal component scores:

```python
_, scores, _ = hl.hwe_normalized_pca(
    for_pca.GT, k=max(1, min(sample_num // 3, 10)), compute_loadings=False
)
```

With big sample sizes, the code above will restrict `k` to the maximum of 10 principal components.

The scores are one of the inputs of the [`hl.pc_relate()` function](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) that will estimate the relatedness between samples in a pairwise manner. To speed things up, I will estimate only the kinship statistic (`statistics="kin"`, while the default is `statistics="all"`. Check the function documentation in the previous link for more details).

```python
relatedness_ht = hl.pc_relate(
    for_pca.GT,
    min_individual_maf=0.01,
    scores_expr=scores[for_pca.col_key].scores,
    block_size=4096,
    min_kinship=0.05,
    statistics="kin",
)
```

The `relatedness_ht` object is a Hail `Table`. It differs from a `MatrixTable` by not having column nor entries fields. We determine related samples by filtering this table to keep only the samples above the kinship threshold and then passing the object to [`hl.maximal_independent_set()` function](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set):

```python
pairs = relatedness_ht.filter(relatedness_ht["kin"] > RELATEDNESS)

related_samples_to_remove = hl.maximal_independent_set(pairs.i, pairs.j, False)
```

The `related_samples_to_remove` object will contain a selection of samples to be removed from the dataset because they come from related individuals in the sample. We perform this filtering with the command below. Notice the use of the keyword `keep=False` to *negate* the selection (I do *not* want to keep related samples).

```python
mt = mt.filter_cols(hl.is_defined(related_samples_to_remove[mt.col_key]), keep=False)
```

Use `mt.count_cols()` to assess if any sample was removed.

### Wrapping up the sample QC

Here I finish the sample QC process. I then save the `relatedness_ht` and the dataset object `mt` to disk with `write()`:

```python
relatedness_ht.write("relatedness.ht", overwrite=True)

mt.write("sampleqc_pass.mt", overwrite=True)
```

I can now proceed to variant QC.

## Variant QC

I init Hail again if needed, import some more functions from gnomAD, create variables with variant quality thresholds, read the `MatrixTable` with the data passing sample QC and calculate common variant statistics (such as allelic frequency) with [`hl.variant_qc()`](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.variant_qc):

```python
# hail_demo_variant_qc.py
# Import modules and init Hail
import hail as hl

from gnomAD.utils.annotations import bi_allelic_site_inbreeding_expr
from gnomAD.variant_qc.random_forest import apply_rf_model, median_impute_features
from gnomAD.variant_qc.pipeline import train_rf_model
from hail_init import DEFAULT_REF

# Variant Quality hard filters
INBR_COEFF = -0.3
AB_LOWER_LIM = 0.2
AB_UPPER_LIM = 1 - AB_LOWER_LIM

# Read MatrixTable with sample QC-passing dataset
mt = hl.read_matrix_table("sampleqc_pass.mt")

mt = hl.variant_qc(mt)
```

Even though allelic frequency may be already been lift over from the original VCF input, I recommend using `hl.variant_qc()` since it calculates potentially useful information besides allelic frequency, such as p-values from the test of Hardy-Weinberg equilibrium. Check the function documentation at Hail to see the complete list of statistics.

### Filter variants by genotype-related quality thresholds

Next, I calculate two variant-level metrics: inbreeding coefficient and the maximum p-value for sampling the observed allele balance under a binomial model.

```python
mt = mt.annotate_rows(inbr_coeff=bi_allelic_site_inbreeding_expr(mt.GT))

mt = mt.annotate_rows(
    pab_max=hl.agg.max(
        hl.binom_test(mt.AD[1], mt.DP, 0.5, "two-sided")
    )
)
```

Next, I remove variants with excess of heterozygotes by inbreeding coefficient and variants for which no sample had high-quality genotypes by evaluating allele balance (the proportion of sequencing reads that support the variant):

```python
mt = mt.filter_rows(mt.inbr_coeff > INBR_COEFF)

# Removing variants for which no sample had high quality genotypes with hl.any()
mt = mt.filter_rows(hl.agg.any(mt.GQ >= 20))
mt = mt.filter_rows(hl.agg.any(mt.DP >= 10))

mt = mt.annotate_entries(AB=(mt.AD[1] / hl.sum(mt.AD))) # AB = allele balance

mt = mt.filter_rows(
    hl.agg.any(
        (mt.GT.is_hom_ref() & (mt.AB < AB_LOWER_LIM))
        | (mt.GT.is_het() & (mt.AB >= AB_LOWER_LIM) & (mt.AB <= AB_UPPER_LIM))
        | (mt.GT.is_hom_var() & (mt.AB > AB_UPPER_LIM))
    )
)
```

### Variant QC by random forest model

The gnomAD team adopted a [random forest model](https://gnomAD.broadinstitute.org/news/2018-10-gnomAD-v2-1/) to filter out sequencing artifacts. Briefly, they labelled variants passing quality thresholds (such as GATK's VQSR) as true positives and variants not passing as false positives. Next, they performed a supervised random forest training with some variant-level features. From now on, I try to replicate their method with the best of my understanding.

I label the variants as true and false positives:

```python
mt = mt.annotate_rows(tp=hl.if_else(hl.len(mt.filters) == 0, True, False))
mt = mt.annotate_rows(fp=hl.if_else(hl.len(mt.filters) != 0, True, False))
```

The logic is that if the `filters` field (which is carried over from the GenomicsDB-exported VCF) is empty, it indicates the variant passed the VQSR filter and is a false positive otherwise. Thus, I create two Boolean-type columns indicating it. Next, I create a Hail `Table` extracting the needed features and the `tp` and `fp` fields and assigning it to the `rf_ht` object:

```python
rf_ht = mt.select_rows(
    mt.inbr_coeff,
    mt.info.SOR,
    mt.info.ReadPosRankSum,
    mt.info.MQRankSum,
    mt.info.QD,
    mt.pab_max,
    mt.variant_type,
    mt.n_alleles,
    mt.mixed_site,
    mt.spanning_deletion,
    mt.tp,
    mt.fp,
).rows()
```

Remember that most of these features are brought over from the VCF INFO field, while the others were generated with the help of Hail. I also generate a Python list with the features names:

```python
features = [
    "inbr_coeff",
    "SOR",
    "ReadPosRankSum",
    "MQRankSum",
    "QD",
    "pab_max",
    "variant_type",
    "n_alleles",
    "mixed_site",
    "spanning_deletion",
]
```

Since random forest models do not tolerate missing data, I use a gnomAD function to impute any missing data with the median of the field:

```python
rf_ht = median_impute_features(rf_ht)
```

I will reserve all variants located in chromosome 20 to perform model evaluation with the help of the [`hl.parse_locus_interval()` function](https://hail.is/docs/0.2/functions/genetics.html#hail.expr.functions.parse_locus_interval):

```python
test_intervals = ["chr20"]

test_intervals = [
    hl.parse_locus_interval(x, reference_genome="GRCh38") for x in test_intervals
]
```

I now may train the model with the help of gnomAD's `train_rf_model()` function. Internally, the function will select a balanced dataset of true positives and false positives to train the model.

```python
rf_trained_ht, rf_model = train_rf_model(
    rf_ht,
    rf_features=features,
    tp_expr=rf_ht.tp,
    fp_expr=rf_ht.fp,
    test_expr=hl.literal(test_intervals).any(
        lambda interval: interval.contains(rf_ht.locus)
    ),
)
```

The `rf_trained_ht` output is a Hail `Table` with annotations related with the random forest model training. The `rf_model` object is the model binary generated by Spark. The inputs include the `rf_ht Table`, the `features` list, the `rf_ht.tp` and `rf_ht.fp` Boolean columns and a `test_expr` argument receives a expression that will ensure that the loci contained in the interval object `test_intervals` will be used for model evaluation.

After model training, I left join the `rf_ht` with the model-annotated `rf_trained_ht` into the `ht Table`. I use it as the input for the gnomAD's `apply_rf_model()` function. It will apply the random forest model in all variants in the genome, including those not selected for training.

```python
ht = rf_ht.join(rf_trained_ht, how="left")

rf_results = apply_rf_model(
    ht=ht,
    rf_model=rf_model,
    features=hl.eval(rf_trained_ht.features),
    label="rf_label",
    prediction_col_name="rf_prediction",
)
```

I then write to disk the Hail `Table` containing a summary of the number of variants originally labeled as true or false positives and the prediction by the model. In other words - a confusion matrix:

```python
rf_summary_ht = rf_results.group_by(
    "tp", "fp", "rf_train", "rf_label", "rf_prediction"
).aggregate(n=hl.agg.count())

rf_summary_ht.write("rf_summary.ht", overwrite=True)
```

I unpack the `rf_results Table` fields and join them in the sample QC `MatrixTable`:

```python
variantqc_pass = mt.annotate_rows(**rf_results[mt.locus, mt.alleles])
```

It can be easily filtered to keep only variants predicted to be true positives by the model and then written to disk:

```python
variantqc_pass = hl.filter_rows(variantqc_pass.rf_prediction == "TP")

variantqc_pass.write("variantqc_pass.mt", overwrite=True)
```

## Filter loci by coordinates or allelic frequency

Now the dataset has passed all QC, I may query it to search for new variants or answer other scientific questions. To illustrate that, I will filter the dataset to contain variants in delimited regions of the genome with a certain range of allelic frequencies. To parse specific regions from a genome, I can create a Python list of strings representing exact or approximate coordinates:

```python
intervals = ["chr10:52765380-52772784", "chr1:100M-200M"]
```

Now I will apply the filter to a `MatrixTable` with `hl.parse_locus_interval()`:

```python
filtered_mt = hl.filter_intervals(
    variantqc_pass,
    [hl.parse_locus_interval(x, reference_genome=DEFAULT_REF) for x in intervals]) 
```

I can also pinpoint an individual locus and create a window of nucleotides before and after it. In the code below I create a window of 100,000 nucleotides before and after a specific position in chromosome X.

```python
locus = hl.parse_locus("chrX:23833353", DEFAULT_REF)
window = locus.window(100000, 100000)

filtered_mt = variantqc_pass.filter_rows(window.contains(variantqc_pass.locus))

```

If I wanted to check the first few genotypes of the filtered `MatrixTable` by the specified window I would use the command below:

```python
filtered_mt.GT.show()
```

Since this coordinate is not on the pseudoautosomal region (PAR) of the X chromosomes, karyotypically normal male individuals will be haploid around this region, and Hail would correctly show only one allele instead of two in the `GT` call.

Since I used `hl.variant_qc()` at the beginning of variant QC, I may filter the variants by their allelic frequency. For example, If I wanted to keep only the variants with less than 1% frequency I would do:

```python
filtered_mt = filtered_mt.filter_rows(filtered_mt.variant_qc.AF[1] < 0.01)
```

Remember to `write()` the `MatrixTable` if you want to save the dataset with any applied filters.

## Conclusion

In this post I:

- Introduced the Hail framework;
- Mentioned software used for preparing multi-sample VCF input starting with multiple gVCF files;
- Demonstrated how to perform sample and variant QC with Hail;
- Demonstrated how filter dataset according to quality metrics, locus or loci interval.

## Acknowledgements

This demonstration was made possible with the help of by insights acquired by reading through the:

- [Hail Discussion Forum](https://discuss.hail.is/) posts;
- [Centre for Population Genomics GitHub repository](https://github.com/populationgenomics/joint-calling);
- [gnomADs' utilities GitHub repository](https://github.com/broadinstitute/gnomad_methods);
- gnomAD team's supplementary material from their [2020 Nature paper](https://www.nature.com/articles/s41586-020-2308-7).

## Appendix

In this post I gave general directions how to combine multiple gVCF files into one single VCF input. Hail actually has an experimental function that has this very purpose: [`hl.experimental.run_combiner()`](https://hail.is/docs/0.2/experimental/vcf_combiner.html). However, I tried to use this function and had problems with it. It generates a "sparse" `MatrixTable` and unfortunately I found the function documentation insufficiently clear on how to work with this slightly different form of intermediate input, so I resorted to GATK's `GenomicsDBImport` as stated. Since Hail is in active development, I expect improvement on both the function and on its documentation.

Throughout the demonstration I used `write()` method to write `MatrixTable`s to disk and later read them back into the session with `read_matrix_table()`. Alternatively I could have used Hail's `checkpoint()` method as an alias for these sequential operations. Read the documentation [here](https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.checkpoint).

*Subscribe to my [RSS feed](https://antoniocampos13.github.io/feeds/all.rss.xml), [Atom feed](https://antoniocampos13.github.io/feeds/all.atom.xml) or [Telegram channel](https://t.me/joinchat/AAAAAEYrNCLK80Fh1w8nAg) to keep you updated whenever I post new content.*

## References

[Hospital Israelita Albert Einstein](https://www.einstein.br/Pages/Home.aspx)

[Hail |  Index](https://hail.is/)

[Broad Institute](https://www.broadinstitute.org/)

[gnomAD](https://gnomAD.broadinstitute.org/)

[Hail |  Index](https://hail.is/#install)

[gnomAD module | PyPi](https://pypi.org/project/gnomAD/)

[Genomic Analysis Toolkit](https://gatk.broadinstitute.org/hc/en-us)

[bcftools](http://samtools.github.io/bcftools/bcftools.html)

[Setting Up Your Unix Computer for Bioinformatics Analysis](https://antoniocampos13.github.io/setting-up-your-unix-computer-for-bioinformatics-analysis.html)

[Genomic Variant Call Format (gVCF)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format)

[Illumina DRAGEN Bio-IT Platform| Variant calling & secondary genomic analysis software tool](https://www.illumina.com/products/by-type/informatics-products/dragen-bio-it-platform.html)

[Variant Call Format specification | version 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)

[GATK | GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)

[GATK | Variant Qualit Score Recalibration (VQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)

[GATK Resource Bundle at Google Cloud](https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false)

[gsutil tool | Google Cloud](https://cloud.google.com/storage/docs/gsutil)

[Hail | Hail 0.2](https://hail.is/docs/0.2/index.html)

[Hail | MatrixTable Tutorial](https://hail.is/docs/0.2/tutorials/07-matrixtable.html)

[Hail | Genetics | sample_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc)

[Bokeh documentation](https://docs.bokeh.org/en/latest/index.html)

[Hail | Genetics | hwe_normalized_pca](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca)

[Hail | Relatedness](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate)

[Hail | Miscellaneous](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set)

[Hail | Genetics | variant_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.variant_qc)

[gnomAD v2.1 | gnomAD news](https://gnomAD.broadinstitute.org/news/2018-10-gnomAD-v2-1/)

[Hail | Genetics functions](https://hail.is/docs/0.2/functions/genetics.html#hail.expr.functions.parse_locus_interval)

[Hail Discussion](https://discuss.hail.is/)

[Centre for Population Genomics GitHub repository](https://github.com/populationgenomics/joint-calling)

[gnomADs' utilities GitHub repository](https://github.com/broadinstitute/gnomad_methods)

[The mutational constraint spectrum quantified from variation in 141,456 humans](https://www.nature.com/articles/s41586-020-2308-7)

[Hail | VCF Combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

[Hail | MatrixTable | checkpoint](https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.checkpoint)
