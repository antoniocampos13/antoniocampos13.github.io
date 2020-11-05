---
Title: Machine Learning with Python: Supervised Classification of TCGA Prostate Cancer Data (Part 1 - Making Features Datasets) 
Status: draft
Date: 2020-11-06 15:00
Author: Antonio Victor Campos Coelho
Categories: Python
Tags: Bioinformatics, gene expression, machine learning, supervised classification
---

## Introduction

In a [previous post](https://antoniocampos13.github.io/working-with-cancer-genomics-cloud-datasets-in-a-postgresql-database-part-1.html), I showed how to retrieve [The Cancer Genome Atlas (TCGA)](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) data from the [Cancer Genomics Cloud (CGC) platform](http://www.cancergenomicscloud.org/). I downloaded gene expression quantification data, created a relational database with PostgreSQL, and created a dataset uniting the raw quantification data for 675 differentially expressed genes [identified by edgeR](https://antoniocampos13.github.io/differential-expression-analysis-with-edger-in-r.html#differential-expression-analysis-with-edger-in-r), race, age at diagnosis and tumor size at diagnosis.

In this post, I will use Python to prepare features datasets to use them to produce a classification model using machine learning tools, especially the `scikit-learn` module. Check its documentation [here](https://scikit-learn.org/stable/).

The dataset and code presented here are available in my [portfolio](https://github.com/antoniocampos13/portfolio/tree/master/Python/2020_11_05_Supervised_Machine_Leaning_TCGA_Data).

## Prepare workspace

Since I will generate a statistical model, it is important to test it in never-seen before data, to assess its prediction validity. Thus, I will use a customized script, `make_features.py`, to transform and split the dataset into separate train and test datasets. I will fit the model and then use it to predict the risk of prostate cancer of the subjects included in the test dataset and then compare with actual status (control: prostate cancer in remission/cases: prostate cancer progressing) to see how well it predicted.

Check the `make_features.py`. Let's examine it. First, I import the necessary modules:

```python
from pathlib import Path

import janitor as jn
import pandas as pd
from sklearn import preprocessing, model_selection
```

`pathlib` is a module from the Python standard library. It helps file path management. `janitor` is a module for data clean-up and manipulation. `pandas` is also used for data manipulation and transformation, especially if it is contained on data frames. `sklearn` is an alias for `scikit-learn`. It is a powerhouse of several machine learning functions and utilities. Install all non-standard library modules using your Python package manager, usually `pip` or through anaconda, if you are using a conda environment with Python.

Next, I set up some constants that I will use later:

```python
RANDOM_SEED=42
TEST_SIZE=0.30
INDEX=676
CORR_THRESHOLD=0.80
```

* `RANDOM_SEED`: An integer used to initialize the pseudo-random number generator. It helps generate reproducible outputs in different sessions/computers;
* `TEST_SIZE`: A decimal (float) number indicating the size of the test dataset. Currently, it is 30% of the samples in the complete dataset;
* `INDEX`: An integer indicating a slicing index. I will explain it later.
* `CORR_THRESHOLD`: A float number indicating a threshold for eliminating correlated variables in the dataset; it helps overfitting by reducing the [multicollinearity issue](https://en.wikipedia.org/wiki/Multicollinearity).

Next, I will indicate the path of the dataset, which i saved as a CSV file:

```python
project_folder = Path().resolve().parent.parent
CSV_FILE = project_folder / "data" / "interim" / "prostate_cancer_dataset.csv"
```

See below the structure of the current working directory:

```text
.
├── data
│   ├── interim
│   │   └── prostate_cancer_dataset.csv
│   └── processed
├── models
└── src
    ├── features
    │   └── make_features.py
    └── models
        └── make_model.py
```

## Load data into Pandas

Now, I create a `pandas.DataFrame` using the `read_csv` function, and assign it to `df` object:

```python
df = pd.read_csv(CSV_FILE)
```

Check the column names with the `columns` method:

```python
df.columms
```

The `status` column contains the sample labels: `0` means "control" and `1` means "case". The other columns are the variables (or features) of each sample.

## Dummy-encode categorical variables

I can now transform the data, making it suitable for machine learning modeling. Luckily, the TCGA dataset has no missing values and no gross errors are present. Otherwise, I would have to impute missing values and correct the errors. Now I will recode the categorical variables in the dataset. To check which variables are categorical, use the `pandas` `dtypes` method. The variables labeled with `object` are usually categorical. If you see `object` beside the name of a variable you know it is not categorical, then it is possible that there are missing data or some other error.

```python
df.dtypes
df = pd.get_dummies(df, drop_first=True)
```

## Split datasets

With the transformed dataset, let's go ahead and split the dataset into training and testing datasets.

```python
X, y = jn.get_features_targets(df, target_columns="status")

X_train, X_test, y_train, y_test = model_selection.train_test_split(
    X, y, test_size=TEST_SIZE, stratify=y, random_state=RANDOM_SEED
)
```

The `get_features_targets` from `janitor` module will get the `df` object and separate the sample labels from the variables (features). The features will then be assigned to object `X`, which will be instanced as a `pandas.DataFrame` and the `status` column (see the `target_columns` argument) will be put in the `y` object, which will be instanced as a `pandas.Series`.

The second function, `model_selection.train_test_split` from `sklearn` will split the `X` and `y` objects into its training and testing counterparts -- therefore, four objects: `X_train`, `X_test`, `y_train` and `y_test`. Check the `TEST_SIZE` and `RANDOM_SEED` constants I set up in the beginning of the script being used here. They indicate that 30% of the dataset must be included as a test dataset (therefore 70% in the training dataset) and setting a integer into the `random_state` argument ensures that the splitting outputs can be reproduced.

Also note the `stratify` argument. It ensures that the same proportion of cases and controls are drawn for training and testing datasets. For this, I indicate the object containing the labels, which is `y`.

The commands below print the number of cases and controls for each dataset:

```python
print(f"Training Set has {sum(y_train)} Positive Labels (cases) and {len(y_train) - sum(y_train)} Negative Labels (controls)")
print(f"Test Set has {sum(y_test)} Positive Labels (cases) and {len(y_test) - sum(y_test)} Negative Labels (controls)")
```

The output:

```text
Training Set has 38 Positive Labels (cases) and 127 Negative Labels (controls)
Test Set has 16 Positive Labels (cases) and 55 Negative Labels (controls)
```

## Normalize data

Now I will standardize (scale or [Z-score normalize](https://en.wikipedia.org/wiki/Standard_score)) the numerical columns. For each numerical feature, I will calculate its mean and standard deviation. Then, for each observed value of the variable, I will subtract the mean and divide by the standard deviation.

After the dummy coding of the categorical variables, they were transposed to the end of the data frame. I manually checked the column names and identified the column index. I then sliced the list of column names and stored in the `num_cols` object:

```python
num_cols = list(X.columns)[:INDEX]
```

That's why I conveniently saved the index number into the `INDEX` variable at the beginning of the script. It helps code reusability with different data. It is one of the advantages of avoiding the so called ["magic numbers"](https://www.pluralsight.com/tech-blog/avoiding-magic-numbers/).

With the numeric columns identified, I then used the `StandardScaler()` function from `sklearn`'s `preprocessing` module:

```python
sca = preprocessing.StandardScaler()
X_train[num_cols] = sca.fit_transform(X_train[num_cols])
X_test[num_cols] = sca.transform(X_test[num_cols])
```

I used the function `fit_transform()` function in the `X_train` dataset, passing the `num_cols` list as the indication of which columns need to be normalized (using the brackets `[]` notation) and saving the resulting transformation with the same names into the `X_train` object. Thus, for all effects and purposes, I am replacing the old columns with the normalized columns.

Note that for `X_test` dataset I used a different formula, `transform()`. This is because I am using just the coefficients fitted by `fit_transform()` in the train dataset to generate the normalization in the test dataset. This way I am sure that the scaling is not "contaminated" by the test data, that is supposed to not seem before the classification model fitting.

## Remove correlated features

Now I will filter out correlated features. Please note that I would not normally do this with genomic data, but since here is just an exercise, I will show how to do it so you can apply to your projects when necessary. Check below the code (hat tip to [Dario Radečić](https://towardsdatascience.com/feature-selection-in-python-recursive-feature-elimination-19f1c39b8d15) for the snippet):

```python
correlated_features = set()
correlation_matrix = X_train.corr()

for i in range(len(correlation_matrix.columns)):
    for j in range(i):
        if abs(correlation_matrix.iloc[i, j]) > CORR_THRESHOLD:
            colname = correlation_matrix.columns[i]
            correlated_features.add(colname)
```

The commands above will calculate the pairwise correlation between all features in the dataset, creating a `correlation_matrix`. If the correlation between two features is above `CORR_THRESHOLD` (currently 0.80), the loop will store the name of one of them into the `correlated_features` set, ensuring that no names are repeated. If you want to know how may pairs of correlated features were present in the dataset, run the command below to print to the console:

```python
print(f"Number of correlated feature pairs: {len(correlated_features)}")
```

Then I can convert the set as list and pass it to `pandas`' `drop()` method:

```python
X_train = X_train.drop(list(correlated_features),axis=1)
X_test = X_test.drop(list(correlated_features),axis=1)
```

I dropped the columns and saved the data frames with the same name for convenience. Note the `axis=1` argument, it tells `drop()` to look the *columns* for the list elements and then remove them.

## Serializing datasets for future use

Finally, the train and test datasets are ready for use! To save them into disk for later use, I will use `pandas`' `to_pickle` method. A "pickled" Python object is *serialized*: it was converted into a series of bytes (a byte stream). This byte stream therefore can be "unpickled" to restore the object exactly as it was when was "pickled". This is useful to backup data, share with others, and so on. Using `pathlib.Path` notation, I saved the objects into the "data/processed/" folder:

```python
X_train.to_pickle(project_folder/"data"/"processed"/"X_train")
X_test.to_pickle(project_folder/"data"/"processed"/"X_test")
y_train.to_pickle(project_folder/"data"/"processed"/"y_train")
y_test.to_pickle(project_folder/"data"/"processed"/"y_test")
```

**Note: pickled objects are Python-specific only -- non-Python programs may not be able to reconstruct pickled Python objects. WARNING: never, NEVER, unpickle data you do not trust. As it says in the Python documentation: ["It is possible to construct malicious pickle data which will execute arbitrary code during unpickling"](https://docs.python.org/3/library/pickle.html).**

Now, my folders are like this:

```text
.
├── data
│   ├── interim
│   │   └── prostate_cancer_dataset.csv
│   └── processed
│       ├── X_train
│       ├── X_test
│       ├── y_train
│       └── y_test
├── models
└── src
    ├── features
    │   └── make_features.py
    └── models
        └── make_model.py
```

## Conclusion

In this post, I showed how to:

* Import a CSV dataset into `pandas`;
* Dummy-encode categorical data;
* Split the dataset into train/test datasets;
* Normalize data (Z-scores);
* Serialize (pickle) the datasets for future use.

Go to the [Part 2](), where I show how to use the datasets to generate a classification model for predicting risk of prostate cancer disease progression with the `make_model.py` script.

*Subscribe to my [RSS feed](https://antoniocampos13.github.io/feeds/all.rss.xml), [Atom feed](https://antoniocampos13.github.io/feeds/all.atom.xml) or [Telegram channel](https://t.me/joinchat/AAAAAEYrNCLK80Fh1w8nAg) to keep you updated whenever I post new content.*

## References

[The Cancer Genome Atlas Program](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga)

[Cancer Genomics Cloud](http://www.cancergenomicscloud.org/)

[scikit-learn: machine learning in Python - scikit-learn 0.23.2 documentation](https://scikit-learn.org/stable/)

[Multicollinearity - Wikipedia](https://en.wikipedia.org/wiki/Multicollinearity)

[Standard score - Wikipedia](https://en.wikipedia.org/wiki/Standard_score)

[Avoiding Magic Numbers](https://www.pluralsight.com/tech-blog/avoiding-magic-numbers/)

[Feature Selection in Python — Recursive Feature Elimination](https://towardsdatascience.com/feature-selection-in-python-recursive-feature-elimination-19f1c39b8d15)

[pickle — Python object serialization | Python 3.9.0 documentation](https://docs.python.org/3/library/pickle.html)
