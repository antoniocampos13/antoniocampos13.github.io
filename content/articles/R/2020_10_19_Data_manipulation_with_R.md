---
Title: Data manipulation with R
Status: draft
Date: 2020-10-19 13:30
Author: Antonio Victor Campos Coelho
Categories: R
Tags: Bioinformatics, gene expression, SQL, PostgreSQL
---

## Introduction

In my [previous post](https://antoniocampos13.github.io/working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1.html#working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1) I demonstrated how to obtain a prostate cancer dataset with genomic information in the form of gene expression quantification and created a local PostgreSQL database to hold the data.

Here, I will use R to connect to the PostgreSQL database, retrieve and then prepare the CGC data to perform a Differential Expression Analysis for sequence count data in the CGC dataset. This demonstration is on a RStudio project running in Windows 10, but the same steps can be followed on an Unix system.

As always, the code demonstrated here is available on my [portfolio on GitHub](https://github.com/antoniocampos13/portfolio/tree/master/R/2020_10_19_Data_manipulation_with_R).

## Setting up .Renviron file

In my working directory, I create a text file named `.Renviron` to store the credentials of the PostgreSQL database with the following information:

```r
userid = "<USER_NAME>"
pwd = "<PASSWORD>"
```

Replace `<USER_NAME>` and `<PASSWORD>` with your credentials used to access the PostgreSQL database. Usually, the default username is `postgres` and the password is defined during PostgreSQL installation.

## Install/Load packages

Then I open a RStudio session and create a project in the folder containing the `.Renviron` file. Now, I need to load the packages I will use today. You can install them using `install.packages()` function if you do not have them installed yet:

```r
# Run only once
install.packages(c("RPostgres", "here", "tidyverse"))
```

Note that to install more than one package at once we must use the concatenate `c()` command to pass the packages names and they must be quoted and separated by commas -- a R *vector*. The package dependencies will be installed as well.

Now I load the packages into the R session:

```r
library(RPostgres)
library(here)
library(tidyverse)
```

The `RPostgres` is needed to connect to the PostgreSQL database I created before. The [`here` package](https://github.com/jennybc/here_here) handles file paths. `tidyverse` is a powerful collection of packages for data manipulation. By using it, it will actually load several other packages, such as `dplyr`, `stringr` and `tidyr`. It is one of the packages I use the most. Check `RPostgres` documentation [here](https://www.rdocumentation.org/packages/RPostgres/versions/1.2.1) and `tidyverse` documentation [here](https://www.tidyverse.org/).

## Set up connection to tcga database

I now set up the connection to the database:

```r
con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "tcga",
  host = "localhost",
  port = 5432,
  user = Sys.getenv("userid"),
  password = Sys.getenv("pwd")
)
```

The `Sys.getenv()` will retrieve the information in the `.Renviron` file. Remember to not share this file so your credentials remain secret. The connection credentials are now stored in the `con` object.

## Retrieve and pivot gene_counts_cases table

With the command below, I retrieve the table containing the case-identified, raw gene counts I created last time:  

```r
cases <- dbGetQuery(con, "SELECT * FROM gene_counts_cases")
```

The `dbGetQuery()` is one of the functions from `DBI` package, a dependency of `RPostgres`. Check `DBI` documentation [here](https://www.rdocumentation.org/packages/DBI/versions/0.5-1).

The information is now stored in the `cases` object, which has three columns: `case_id`, `gene_id` and `gene_counts`. Therefore, each row is a combination of a case, a gene and a gene count -- it is a "long" format.

I now will reformat the table to a "wide" format (pivot) because it is a requirement for the differential expression analysis. The pivoted table will have *G*x*N* dimensions, where *G* are the number of rows (the number of gene transcripts quantified) and *N* the number of columns (the number of cases).

I will use the `pivot_wider()` function of the `tidyr` package to do the job (**WARNING: COMPUTATION INTENSIVE STEP**):

```r
cases_pivoted <-
  cases %>% pivot_wider(
    names_from = case_id,
    values_from = gene_count,
    values_fill = 0,
    names_repair = "check_unique",
    values_fn = mean
  )
```

`cases_pivoted` is the name of the object that will hold the pivoted table. The `%>%` is `dplyr`'s syntax. It means that we are piping the contents of `cases` object into the `pivot_wider()` function and its arguments.

The `names_from` argument tells which column will be pivoted to generate new columns. The `values_from` argument tells which column hold the values that will fill the new pivoted table. The `values_fill` argument will substitute any missing data for a zero. The `names_repair` argument checks that each new column has a unique name. Finally, the `values_fn` argument indicates the function that must be applied to the values filling the new pivoted table. Note that I used the `mean` function because during this step I noticed that some cases were associated with more than one gene expression quantification file. Therefore, I had to take the mean of these extra gene counts to correctly generate the pivoted table.

Since I calculated means for the values, I then rounded to the next integer all numerical data in the table with the help of `dplyr`'s `mutate()`, `across()` and `where()` functions and `round()`, which is one of R's standard (base) functions:

```r
counts <- cases_pivoted %>% mutate(across(where(is.numeric), round, 0))
```

The `where(is.numeric)` ensures that I only manipulated numeric data in the table. The `across()` function applies the same transformation (rounding in this case) to multiple columns. Finally, `mutate()` adds the new variables (rounded columns), replacing the existing ones.

## Retrieve and de-duplicate follow_up table

Now I set aside the `counts` table for a moment to prepare the sample classifications. I realized that the `follow_up` table had duplicate data for some reason. Thus I connected to the database and de-duplicate the `followup_primarytherapyoutcomesuccess_1` column by using string aggregation (PostgreSQL's [`STRING_AGG` function](https://www.postgresqltutorial.com/postgresql-aggregate-functions/postgresql-string_agg-function/)). I also changed its name to `outcome` for simplicity, and saved the results as the `followup_dedup` table. This is the command I used:

```r
dbExecute(
  con,
  "CREATE TABLE followup_dedup AS SELECT case_id, STRING_AGG(followup_primarytherapyoutcomesuccess_1, ',') AS outcome FROM follow_up GROUP BY case_id"
)
```

Note how I used `DBI`'s `dbExecute()` instead of `dbGetQuery()`, since I will not retrieve the table into the R session.

Now I create other table to link the cases ID numbers with the de-duplicated outcomes, saving the result as the `outcomes` table:

```r
dbExecute(
  con,
  "CREATE TABLE outcomes AS SELECT case_id, outcome FROM allcases INNER JOIN followup_dedup USING(case_id)"
)
```

I now retrieve the `outcomes` table into the R session:

```r
outcomes <- dbGetQuery(con, "SELECT * FROM outcomes")
```

## Create sample classification and new labels

I will now create a new column named `class` in the outcomes table with a simplified case/control classification based on the `outcome` column with the help of `dplyr`'s `mutate()` and `case_when()`, which vectorizes multiple if/else statements (it is an R equivalent of the SQL `CASE WHEN` statement):

```r
outcomes <- outcomes %>% mutate(class = case_when(
  str_detect(outcome, "Complete") ~ "control",
  str_detect(outcome, "Partial") ~ "case",
  str_detect(outcome, "Disease") ~ "case"
))
```

Now I will create a new label for the cases for simplification. I created a column named `new_names` in the `outcomes` table. The new names were created by joining the classification created in the previous step with the 12 last characters of the `case_id`:

```r
outcomes <- outcomes %>% mutate(new_names = paste0(str_sub(case_id, -12), "_", class))
```

Finally, let's apply the new case labels, substituting the old ones with the help of `dplyr`'s `recode()`  function:

```r
colnames(counts) <- dplyr::recode(colnames(counts), !!!setNames(as.character(outcomes$new_names), outcomes$case_id))
```

## Convert gene_ids into row names, then delete gene_id column

The table is almost in the state required for differential expression analysis. I just need to convert the `gene_id` column into row names of the data frame:

```r
counts <- as.data.frame(counts)

row.names(counts) <- counts$gene_id
```

The `gene_id` is not needed anymore. I can delete it:

```r
counts <- subset(counts, select = -c(gene_id))
```

The data frame is now ready for differential expression analysis for sequence count data. The features (gene IDs) are embedded on the R object row names, whereas each column corresponds to a individual sample. The row/column intersection are therefore, the raw counts of gene expression.

Check the dimensions (*G*x*N*) of the data frame with the `dim()` function and see that there 60483 rows and 236 columns, meaning that 60483 transcripts where quantified in samples obtained from 236 individuals with prostate cancer.

```r
dim(counts)

# Output:
[1] 60483   236
```

See below first few rows and columns of the finalized `counts` data frame:

![Final gene counts table]({static}/images/counts_final.PNG)

In a future post I will demonstrate the differential expression analysis per se.

## Conclusion

* I demonstrated how to connect to a local PostgreSQL database with `Rpostgres` and `DBI` packages;
* Reorganized data by pivoting a long data frame to a wider data frame with gene IDs in rows and samples in columns with `dplyr` package functions;
* Labeled columns as cases or controls for the differential expression analysis, also with `dplyr` package.

*Go back to [Working with Cancer Genomics Cloud datasets in a PostgreSQL database Part 1](https://antoniocampos13.github.io/working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1.html#working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1) and [Part 2](https://antoniocampos13.github.io/working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-2.html#working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-2).*

## References

[Ode to the here package](https://github.com/jennybc/here_here)

[RPostgres package | R Documentation](https://www.rdocumentation.org/packages/RPostgres/versions/1.2.1)

[Tidyverse](https://www.tidyverse.org/)

[DBI package | R Documentation](https://www.rdocumentation.org/packages/DBI/versions/0.5-1)

[PostgreSQL STRING_AGG() Function By Practical Examples](https://www.postgresqltutorial.com/postgresql-aggregate-functions/postgresql-string_agg-function/)
