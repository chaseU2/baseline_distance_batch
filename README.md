# baselinebeta

`baselinebeta` is a small Python package for longitudinal beta-diversity analyses based on **distance to each subject's own baseline sample**.

It can:

- read QIIME 2 distance matrices (`.qza`) without requiring QIIME 2 at runtime;
- read TSV/CSV distance matrices;
- calculate distance-to-individual-baseline tables for one or many beta-diversity metrics;
- define arbitrary metadata subsets;
- automatically generate all subset combinations from chosen metadata columns;
- run Kruskal–Wallis tests per timepoint;
- run pairwise two-sided Mann–Whitney U tests per timepoint;
- report mean/median differences and Cliff's delta;
- calculate Benjamini–Hochberg FDR both across all pairwise tests in a metric/test and separately within each timepoint;
- create multi-panel boxplots with individual subject points and significance brackets;
- run from Python or YAML/CLI.

## Installation

```bash
pip install git+https://github.com/YOUR_USERNAME/baselinebeta.git
```

For development:

```bash
git clone https://github.com/YOUR_USERNAME/baselinebeta.git
cd baselinebeta
pip install -e .
```

## Quick start

```python
from baselinebeta import BaselineDistanceAnalysis

analysis = BaselineDistanceAnalysis(
    metadata="metadata.tsv",
    distances={
        "bray_curtis": "bray-curtis-distance-matrix.qza",
        "unweighted_unifrac": "unweighted-unifrac-distance-matrix.qza",
    },
    sample_col="sample_name",
    subject_col="host_subject_id",
    time_col="timepoint_numeric",
    baseline_time=-1,
    output_dir="baselinebeta_results",
    metric_titles={
        "bray_curtis": "Bray–Curtis",
        "unweighted_unifrac": "Unweighted UniFrac",
    },
    time_labels={
        -1: "Baseline",
        1: "Day 1",
        3: "Day 3",
        7: "Day 7",
    },
)

analysis.add_test(
    name="treatment_all",
    group_col="description_of_treatment",
    levels=["cold", "cold_radiation", "warm", "warm_radiation"],
    comparisons=[
        ("cold", "cold_radiation"),
        ("cold", "warm"),
        ("cold_radiation", "warm_radiation"),
        ("warm", "warm_radiation"),
    ],
)

distance_table = analysis.build_distance_table(
    time_values=[-1, 1, 3, 7],
)

results = analysis.run()

analysis.plot(
    results=results,
    time_order=[1, 3, 7],
    group_colors={
        "cold": "#4C78A8",
        "cold_radiation": "#9ECAE1",
        "warm": "#E45756",
        "warm_radiation": "#F3A6A5",
    },
    group_labels={
        "cold": "Cold",
        "cold_radiation": "Cold + radiation",
        "warm": "Warm",
        "warm_radiation": "Warm + radiation",
    },
)
```

## Automatic subset combinations

The package can generate tests systematically from metadata columns.

```python
analysis.add_tests_from_columns(
    test_cols=["treatment", "mice_model"],
    subset_cols=["treatment", "mice_model", "sex"],
    include_unfiltered=True,
    max_subset_columns=2,
    level_orders={
        "treatment": ["sham", "irradiation"],
        "mice_model": ["WT", "APC"],
    },
)
```

For `test_cols=["treatment"]`, examples generated automatically are:

```text
all samples
mice_model=WT
mice_model=APC
sex=male
sex=female
mice_model=WT & sex=male
mice_model=WT & sex=female
...
```

The tested column itself is automatically excluded from subsetting.

## Output files

A typical output directory contains:

```text
distance_to_baseline_table.tsv
distance_to_baseline_table.csv
kruskal_wallis.tsv
pairwise_mannwhitney.tsv
skipped_baselines.tsv
skipped_tests.tsv
plots/
```

## Statistical interpretation

For each subject, each later sample is compared with that subject's own baseline sample using the chosen beta-diversity distance matrix. At each timepoint, the resulting subject-level distances are compared between metadata groups.

Kruskal–Wallis provides an omnibus comparison across the available groups. Pairwise comparisons use two-sided Mann–Whitney U tests. Cliff's delta is reported in the direction **group 2 minus group 1**.

By default, the baseline timepoint is included in the distance table (distance = 0) but excluded from inferential tests.

If multiple samples exist for the same subject and timepoint, their distances are averaged before statistical testing so that the subject remains the independent observational unit.
