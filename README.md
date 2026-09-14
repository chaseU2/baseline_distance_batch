# baselinebeta

`baselinebeta` is a Python package for longitudinal beta-diversity analyses based on the **distance of each subject to its own baseline sample**.

The package can calculate distance-to-baseline tables, perform statistical comparisons between metadata groups, and create publication-style boxplots with individual subject values and significance annotations.

Repository:

```text
https://github.com/chaseU2/baseline_distance_batch
```

Python import:

```python
from baselinebeta import BaselineDistanceAnalysis
```

Current version:

```text
0.1.0
```

---

## Overview

For every subject, one timepoint is defined as the individual baseline.

For example:

```text
Mouse A

Day -1  <- baseline
Day 1
Day 3
Day 14
Day 30
```

`baselinebeta` calculates the distance between every later sample and the baseline sample from the **same subject**:

```text
distance(Mouse A Day -1, Mouse A Day 1)
distance(Mouse A Day -1, Mouse A Day 3)
distance(Mouse A Day -1, Mouse A Day 14)
distance(Mouse A Day -1, Mouse A Day 30)
```

This is repeated independently for every subject and every supplied beta-diversity distance matrix.

The resulting values can then be compared between metadata groups such as:

```text
sham vs irradiation
WT vs APC
male vs female
```

or within arbitrary metadata subsets such as:

```text
WT only
female only
WT + female
irradiated + APC
```

---

## Features

`baselinebeta` supports:

- longitudinal distance-to-baseline analysis
- individual baseline samples for every subject
- multiple beta-diversity metrics in one analysis
- QIIME 2 `.qza` distance matrices
- TSV, TXT and CSV distance matrices
- arbitrary metadata column names
- arbitrary baseline timepoints
- arbitrary metadata subsets
- multiple simultaneous subset conditions
- automatic generation of metadata subset combinations
- Kruskal-Wallis tests
- pairwise two-sided Mann-Whitney U tests
- Cliff's delta
- mean and median group differences
- Benjamini-Hochberg FDR correction
- multi-panel distance-to-baseline plots
- boxplots with individual subject points
- significance brackets
- configurable colors and labels
- Python API
- YAML configuration
- command-line interface

QIIME 2 itself is **not required at runtime** to read `.qza` distance matrices.

---

# Installation

## Install directly from GitHub

```bash
python -m pip install \
  git+https://github.com/chaseU2/baseline_distance_batch.git
```

Check the installation:

```bash
python -c "import baselinebeta; print(baselinebeta.__version__)"
```

Expected output:

```text
0.1.0
```

---

## Installation on the JLab server

The JLab server may not allow pip to download temporary build dependencies from PyPI.

Activate the existing environment:

```bash
conda activate rachis-qiime2-2026.4
```

Then install without build isolation:

```bash
python -m pip install \
  --no-build-isolation \
  git+https://github.com/chaseU2/baseline_distance_batch.git
```

To force-install the newest version after changes have been pushed to GitHub:

```bash
python -m pip install \
  --upgrade \
  --no-deps \
  --no-build-isolation \
  --force-reinstall \
  git+https://github.com/chaseU2/baseline_distance_batch.git
```

Check the installed version:

```bash
python -c "import baselinebeta; print(baselinebeta.__version__)"
```

---

# Development installation

Clone the repository:

```bash
git clone git@github.com:chaseU2/baseline_distance_batch.git
cd baseline_distance_batch
```

Create a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
```

Install the package in editable mode:

```bash
python -m pip install -e .
```

Run the tests:

```bash
pytest
```

---

# Quick start

```python
from pathlib import Path

from baselinebeta import BaselineDistanceAnalysis


BASE = Path("/path/to/project")


analysis = BaselineDistanceAnalysis(

    metadata=BASE / "metadata.tsv",

    distances={

        "bray_curtis":
            BASE / "bray-curtis-distance-matrix.qza",

        "weighted_unifrac":
            BASE / "weighted-unifrac-distance-matrix.qza",

        "unweighted_unifrac":
            BASE / "unweighted-unifrac-distance-matrix.qza",
    },

    sample_col="sample_name",

    subject_col="host_subject_id",

    time_col="timepoint_numeric",

    baseline_time=-1,

    output_dir=BASE / "baselinebeta_results",

    metric_titles={
        "bray_curtis": "Bray-Curtis",
        "weighted_unifrac": "Weighted UniFrac",
        "unweighted_unifrac": "Unweighted UniFrac",
    },

    time_labels={
        -1: "Baseline",
        1: "Day 1",
        3: "Day 3",
        14: "Day 14",
        30: "Day 30",
    },

    duplicate_baseline="error",
    missing_baseline="skip",
)
```

---

# Build the distance-to-baseline table

```python
distance_table = analysis.build_distance_table(

    time_values=[
        -1,
        1,
        3,
        14,
        30,
    ],

    include_baseline=True,

    save=True,

    verbose=True,
)
```

Each row describes the distance between one sample and the baseline sample of the same subject.

Important output columns include:

```text
metric
subject
sample_id
baseline_sample
timepoint_value
timepoint_label
dist_to_baseline
```

Original metadata columns are retained as well. This allows later subsetting by variables such as:

```text
treatment
mice_model
sex
diet
batch
cage
```

---

# Baseline handling

By default, each subject should have exactly one baseline sample.

Example:

```python
analysis = BaselineDistanceAnalysis(
    ...,
    baseline_time=-1,
    duplicate_baseline="error",
    missing_baseline="skip",
)
```

## Duplicate baseline samples

```text
duplicate_baseline="error"
```

raises an error if a subject has more than one baseline sample.

Alternatively:

```text
duplicate_baseline="first"
```

uses the first available baseline sample.

## Missing baseline samples

```text
missing_baseline="skip"
```

skips subjects without a baseline sample.

Alternatively:

```text
missing_baseline="error"
```

raises an error.

Skipped subjects are documented in:

```text
skipped_baselines.tsv
```

---

# Define statistical tests manually

A test compares distance-to-baseline values between metadata groups at each timepoint.

## Example: treatment across all subjects

```python
analysis.add_test(

    name="treatment_all",

    group_col="treatment",

    levels=[
        "sham",
        "irradiation",
    ],

    comparisons=[
        ("sham", "irradiation"),
    ],

    run_kruskal=False,
    run_pairwise=True,
)
```

---

## Example: treatment within WT subjects

```python
analysis.add_test(

    name="treatment_within_WT",

    group_col="treatment",

    subset={
        "mice_model": "WT",
    },

    levels=[
        "sham",
        "irradiation",
    ],

    comparisons=[
        ("sham", "irradiation"),
    ],
)
```

This means:

```text
first subset:
mice_model == WT

then test:
sham vs irradiation
```

---

## Multiple subset conditions

Multiple metadata filters can be combined:

```python
analysis.add_test(

    name="treatment_within_WT_female",

    group_col="treatment",

    subset={
        "mice_model": "WT",
        "sex": "female",
    },

    levels=[
        "sham",
        "irradiation",
    ],

    comparisons=[
        ("sham", "irradiation"),
    ],
)
```

This means:

```text
mice_model == WT
AND
sex == female
```

before the treatment comparison is performed.

---

# Automatically generate subset combinations

Tests can also be generated automatically from metadata columns.

```python
analysis.add_tests_from_columns(

    test_cols=[
        "treatment",
        "mice_model",
    ],

    subset_cols=[
        "treatment",
        "mice_model",
        "sex",
    ],

    include_unfiltered=True,

    max_subset_columns=2,

    level_orders={
        "treatment": [
            "sham",
            "irradiation",
        ],

        "mice_model": [
            "WT",
            "APC",
        ],
    },
)
```

For a treatment test this can automatically generate:

```text
treatment | all subjects

treatment | mice_model=WT
treatment | mice_model=APC

treatment | sex=male
treatment | sex=female

treatment | mice_model=WT & sex=male
treatment | mice_model=WT & sex=female

treatment | mice_model=APC & sex=male
treatment | mice_model=APC & sex=female
```

The column currently being tested is automatically excluded from subset generation.

For example:

```text
test column:
treatment
```

will not generate:

```text
subset:
treatment=sham
```

because only one treatment level would remain.

---

# Pairwise comparisons

## All possible pairwise comparisons

```python
analysis.add_test(

    name="treatment",

    group_col="treatment",

    levels=[
        "cold",
        "cold_radiation",
        "warm",
        "warm_radiation",
    ],

    comparisons="all",
)
```

This automatically tests all pairwise combinations.

---

## Selected pairwise comparisons

You can also define only biologically relevant comparisons:

```python
analysis.add_test(

    name="treatment",

    group_col="treatment",

    levels=[
        "cold",
        "cold_radiation",
        "warm",
        "warm_radiation",
    ],

    comparisons=[
        ("cold", "cold_radiation"),
        ("cold", "warm"),
        ("cold_radiation", "warm_radiation"),
        ("warm", "warm_radiation"),
    ],
)
```

---

# Run the statistical analysis

```python
results = analysis.run(

    distance_table=distance_table,

    include_baseline_tests=False,

    save=True,

    verbose=True,
)
```

The returned result object contains:

```python
results.distance_table
results.kruskal
results.pairwise
results.skipped
```

For example:

```python
display(results.pairwise)
```

---

# Kruskal-Wallis tests

Kruskal-Wallis can be used as an omnibus test when more than two groups are present.

The null hypothesis is that the distance-to-baseline distributions do not differ between the tested groups.

Example:

```python
analysis.add_test(

    name="treatment",

    group_col="treatment",

    levels=[
        "cold",
        "cold_radiation",
        "warm",
        "warm_radiation",
    ],

    run_kruskal=True,
)
```

Kruskal-Wallis can also be disabled:

```python
analysis.add_test(
    name="treatment",
    group_col="treatment",
    run_kruskal=False,
)
```

---

# Mann-Whitney U tests

Pairwise comparisons use two-sided Mann-Whitney U tests.

Important output columns include:

```text
metric
test_name
timepoint_value
timepoint_label

group_1
group_2

n_1
n_2

mean_1
mean_2
mean_difference_2_minus_1

median_1
median_2
median_difference_2_minus_1

percent_median_difference_2_minus_1

cliffs_delta_2_minus_1

U_statistic
p_value

q_value_fdr_bh_global
q_value_fdr_bh_per_timepoint
```

---

# Cliff's delta

Cliff's delta is reported in the direction:

```text
group 2 minus group 1
```

Therefore:

```text
positive Cliff's delta
```

means that group 2 tends to have larger distance-to-baseline values.

```text
negative Cliff's delta
```

means that group 1 tends to have larger distance-to-baseline values.

The possible range is:

```text
-1 to +1
```

---

# Mean and median differences

The following columns are reported:

```text
mean_difference_2_minus_1
median_difference_2_minus_1
percent_median_difference_2_minus_1
```

These are calculated as:

```text
group 2 - group 1
```

For example:

```text
group_1 = WT
group_2 = APC
```

A positive median difference means:

```text
APC has a larger median distance to baseline than WT
```

---

# FDR correction

Pairwise Mann-Whitney U tests include two Benjamini-Hochberg corrections.

## Global within metric and test

```text
q_value_fdr_bh_global
```

This corrects across the pairwise tests and timepoints belonging to the same:

```text
metric × test definition
```

For example:

```text
Bray-Curtis
×
treatment_within_WT
```

---

## Per timepoint

```text
q_value_fdr_bh_per_timepoint
```

This corrects across pairwise comparisons performed at one timepoint within the same metric and test definition.

If only one comparison exists at a timepoint, the per-timepoint q-value is identical to the raw p-value.

---

# Plotting

`baselinebeta` can generate multi-panel distance-to-baseline plots.

Example:

```python
plot_files = analysis.plot(

    results=results,

    output_dir="baselinebeta_results/plots",

    time_order=[
        1,
        3,
        14,
        30,
    ],

    time_labels={
        1: "Day 1",
        3: "Day 3",
        14: "Day 14",
        30: "Day 30",
    },

    q_column=
        "q_value_fdr_bh_global",

    alpha=0.05,

    show_only_significant=True,

    group_colors={
        "WT": "#DD8452",
        "APC": "#4C72B0",
    },

    group_labels={
        "WT": "WT",
        "APC": "APC",
        "sham": "Sham",
        "irradiation": "Irradiation",
    },

    save_png=True,
    save_pdf=True,

    show=True,
)
```

Each panel represents one timepoint.

The plots contain:

- boxplots
- individual subject values
- black borders around individual points
- sample-size labels
- pairwise significance brackets
- FDR-adjusted significance stars

The y-axis represents:

```text
Distance to individual baseline
```

---

# Significance annotations

By default, significance stars follow:

```text
q < 0.05   *
q < 0.01   **
q < 0.001  ***
```

If:

```python
show_only_significant=True
```

only FDR-significant comparisons are shown.

If:

```python
show_only_significant=False
```

non-significant comparisons can also be displayed.

---

# Statistical interpretation

Suppose Day -1 is defined as the individual baseline.

For subject A at Day 14:

```text
distance_to_baseline
=
distance(
    subject A at Day -1,
    subject A at Day 14
)
```

Now assume the groups are:

```text
WT
vs
APC
```

The Mann-Whitney U test asks:

> Do WT and APC subjects differ in how far their microbiome has moved away from their own baseline by Day 14?

This is different from PERMANOVA.

PERMANOVA asks whether the overall multivariate microbiome composition differs between groups.

`baselinebeta` instead asks whether groups differ in the **magnitude of within-subject change from their own baseline**.

---

# Repeated samples at the same timepoint

If one subject has multiple samples at the same timepoint, their distance-to-baseline values are averaged before statistical testing.

Therefore one subject contributes one independent value per timepoint.

---

# Supported distance-matrix formats

The following inputs are supported:

```text
.qza
.tsv
.txt
.csv
pandas.DataFrame
```

A QIIME 2 installation is not required to read `.qza` distance matrices.

---

# Example: three beta-diversity metrics

```python
analysis = BaselineDistanceAnalysis(

    metadata="metadata.tsv",

    distances={

        "bray_curtis":
            "bray-curtis-distance-matrix.qza",

        "weighted_unifrac":
            "weighted-unifrac-distance-matrix.qza",

        "unweighted_unifrac":
            "unweighted-unifrac-distance-matrix.qza",
    },

    sample_col="sample_name",

    subject_col="host_subject_id",

    time_col="timepoint_numeric",

    baseline_time=-1,
)
```

All analyses are then performed separately for each beta-diversity metric.

---

# Typical output files

A typical output directory contains:

```text
baselinebeta_results/
├── distance_to_baseline_table.tsv
├── distance_to_baseline_table.csv
├── kruskal_wallis.tsv
├── pairwise_mannwhitney.tsv
├── skipped_baselines.tsv
├── skipped_tests.tsv
└── plots/
    ├── ...
    └── ...
```

---

# Command-line usage

`baselinebeta` can also be run with a YAML configuration file.

Example:

```bash
baselinebeta analysis.yaml
```

Run without plots:

```bash
baselinebeta analysis.yaml --no-plots
```

An example configuration is available under:

```text
examples/coldshield_example.yaml
```

---

# Updating the package

After pushing changes to GitHub, reinstall the newest version on the JLab server:

```bash
conda activate rachis-qiime2-2026.4

python -m pip install \
  --upgrade \
  --no-deps \
  --no-build-isolation \
  --force-reinstall \
  git+https://github.com/chaseU2/baseline_distance_batch.git
```

Then verify:

```bash
python -c "import baselinebeta; print(baselinebeta.__version__)"
```

---

# Development workflow

After modifying the package locally:

```bash
git add .
git commit -m "Describe changes"
git push
```

Repository:

```text
git@github.com:chaseU2/baseline_distance_batch.git
```

---

# License

MIT
