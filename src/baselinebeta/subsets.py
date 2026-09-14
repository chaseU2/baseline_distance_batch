from __future__ import annotations

from itertools import product
from typing import Any, Iterable

import pandas as pd


def apply_subset(df: pd.DataFrame, subset: dict[str, Any] | None) -> pd.DataFrame:
    out = df.copy()
    if not subset:
        return out

    for column, value in subset.items():
        if column not in out.columns:
            raise KeyError(f"Subset column '{column}' is not present")

        if isinstance(value, (list, tuple, set)):
            wanted = {str(v) for v in value}
            out = out[out[column].astype(str).isin(wanted)].copy()
        else:
            out = out[out[column].astype(str) == str(value)].copy()

    return out


def subset_description(subset: dict[str, Any] | None) -> str:
    if not subset:
        return "all"

    parts = []
    for column, value in subset.items():
        if isinstance(value, (list, tuple, set)):
            values = ",".join(map(str, value))
            parts.append(f"{column}=[{values}]")
        else:
            parts.append(f"{column}={value}")
    return " & ".join(parts)


def generate_subsets(
    metadata: pd.DataFrame,
    subset_cols: Iterable[str],
    *,
    exclude_cols: Iterable[str] | None = None,
    include_unfiltered: bool = True,
    max_subset_columns: int | None = None,
) -> list[dict[str, Any]]:
    exclude = set(exclude_cols or [])
    columns = [c for c in subset_cols if c not in exclude]

    for column in columns:
        if column not in metadata.columns:
            raise KeyError(f"Subset column '{column}' is not present in metadata")

    options: list[list[Any]] = []
    for column in columns:
        levels = list(pd.unique(metadata[column].dropna().astype(str)))
        try:
            levels = sorted(levels)
        except TypeError:
            pass
        options.append([None] + levels)

    if not columns:
        return [{}] if include_unfiltered else []

    subsets: list[dict[str, Any]] = []
    for values in product(*options):
        subset = {
            col: value
            for col, value in zip(columns, values)
            if value is not None
        }
        if not subset and not include_unfiltered:
            continue
        if max_subset_columns is not None and len(subset) > max_subset_columns:
            continue
        subsets.append(subset)

    return subsets
