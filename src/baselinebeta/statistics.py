from __future__ import annotations

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests


def cliffs_delta(x, y) -> float:
    """Cliff's delta in the direction group 2 (y) minus group 1 (x)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) == 0 or len(y) == 0:
        return np.nan

    greater = 0
    less = 0
    for xi in x:
        greater += np.sum(y > xi)
        less += np.sum(y < xi)
    return float((greater - less) / (len(x) * len(y)))


def _bh(values: pd.Series) -> pd.Series:
    out = pd.Series(np.nan, index=values.index, dtype=float)
    mask = values.notna()
    if mask.any():
        out.loc[mask] = multipletests(
            values.loc[mask].astype(float),
            method="fdr_bh",
        )[1]
    return out


def add_fdr_columns(
    kruskal_df: pd.DataFrame,
    pairwise_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    kw = kruskal_df.copy()
    pw = pairwise_df.copy()

    if not kw.empty:
        kw["q_value_fdr_bh"] = np.nan
        for _, idx in kw.groupby(["metric", "test_name"], dropna=False).groups.items():
            kw.loc[idx, "q_value_fdr_bh"] = _bh(kw.loc[idx, "p_value"])
        kw["significant_fdr"] = kw["q_value_fdr_bh"] < 0.05

    if not pw.empty:
        pw["q_value_fdr_bh_global"] = np.nan
        for _, idx in pw.groupby(["metric", "test_name"], dropna=False).groups.items():
            pw.loc[idx, "q_value_fdr_bh_global"] = _bh(pw.loc[idx, "p_value"])

        pw["q_value_fdr_bh_per_timepoint"] = np.nan
        for _, idx in pw.groupby(
            ["metric", "test_name", "timepoint_value"],
            dropna=False,
        ).groups.items():
            pw.loc[idx, "q_value_fdr_bh_per_timepoint"] = _bh(pw.loc[idx, "p_value"])

        pw["significant_fdr_global"] = pw["q_value_fdr_bh_global"] < 0.05
        pw["significant_fdr_per_timepoint"] = pw["q_value_fdr_bh_per_timepoint"] < 0.05

    return kw, pw


def pair_order(levels: list[str], comparisons: list[tuple[str, str]] | None) -> list[tuple[str, str]]:
    if comparisons is not None:
        return [(str(a), str(b)) for a, b in comparisons]
    return list(combinations([str(x) for x in levels], 2))


def summarize_two_groups(x, y) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mean_1 = float(np.mean(x))
    mean_2 = float(np.mean(y))
    median_1 = float(np.median(x))
    median_2 = float(np.median(y))
    mean_diff = mean_2 - mean_1
    median_diff = median_2 - median_1
    pct_median_diff = (
        median_diff / median_1 * 100.0
        if median_1 != 0
        else np.nan
    )

    return {
        "mean_1": mean_1,
        "mean_2": mean_2,
        "mean_difference_2_minus_1": mean_diff,
        "median_1": median_1,
        "median_2": median_2,
        "median_difference_2_minus_1": median_diff,
        "percent_median_difference_2_minus_1": pct_median_diff,
        "cliffs_delta_2_minus_1": cliffs_delta(x, y),
    }
