from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu

from .distances import load_distance_matrix
from .statistics import add_fdr_columns, pair_order, summarize_two_groups
from .subsets import apply_subset, generate_subsets, subset_description


@dataclass
class BaselineResults:
    distance_table: pd.DataFrame
    kruskal: pd.DataFrame
    pairwise: pd.DataFrame
    skipped: pd.DataFrame


class BaselineDistanceAnalysis:
    """Distance-to-individual-baseline analysis for longitudinal beta diversity."""

    def __init__(
        self,
        *,
        metadata: str | Path | pd.DataFrame,
        distances: dict[str, str | Path | pd.DataFrame],
        sample_col: str,
        subject_col: str,
        time_col: str,
        baseline_time: Any,
        output_dir: str | Path | None = None,
        metric_titles: dict[str, str] | None = None,
        time_labels: dict[Any, str] | None = None,
        duplicate_baseline: str = "error",
        missing_baseline: str = "skip",
    ) -> None:
        self.metadata = self._load_metadata(metadata)
        self.sample_col = sample_col
        self.subject_col = subject_col
        self.time_col = time_col
        self.baseline_time = baseline_time
        self.output_dir = Path(output_dir) if output_dir is not None else None
        self.metric_titles = metric_titles or {}
        self.time_labels = {str(k): str(v) for k, v in (time_labels or {}).items()}
        self.duplicate_baseline = duplicate_baseline
        self.missing_baseline = missing_baseline

        required = [sample_col, subject_col, time_col]
        missing = [c for c in required if c not in self.metadata.columns]
        if missing:
            raise KeyError(f"Metadata is missing required columns: {missing}")

        self.metadata[sample_col] = self.metadata[sample_col].astype("string").str.strip()
        self.metadata[subject_col] = self.metadata[subject_col].astype("string").str.strip()
        self.metadata[time_col] = self.metadata[time_col].astype("string").str.strip()

        if self.metadata[sample_col].duplicated().any():
            dup = self.metadata.loc[
                self.metadata[sample_col].duplicated(False), sample_col
            ].tolist()
            raise ValueError(f"Duplicate sample IDs in metadata: {dup[:20]}")

        if duplicate_baseline not in {"error", "first"}:
            raise ValueError("duplicate_baseline must be 'error' or 'first'")
        if missing_baseline not in {"skip", "error"}:
            raise ValueError("missing_baseline must be 'skip' or 'error'")

        self.distances = {
            name: load_distance_matrix(source)
            for name, source in distances.items()
        }
        self.tests: list[dict[str, Any]] = []

    @staticmethod
    def _load_metadata(source: str | Path | pd.DataFrame) -> pd.DataFrame:
        if isinstance(source, pd.DataFrame):
            return source.copy()

        path = Path(source)
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        md = pd.read_csv(path, sep=sep, dtype=str)
        first_col = md.columns[0]
        md = md[~md[first_col].astype(str).str.startswith("#q2:")].copy()
        return md

    @staticmethod
    def _time_series(series: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().sum() == series.notna().sum():
            return numeric
        return series.astype(str)

    def _time_mask(self, series: pd.Series, value: Any) -> pd.Series:
        numeric = pd.to_numeric(series, errors="coerce")
        try:
            target = float(value)
            if numeric.notna().any():
                return numeric == target
        except (TypeError, ValueError):
            pass
        return series.astype(str) == str(value)

    def _time_label(self, value: Any) -> str:
        key = str(value)
        if key in self.time_labels:
            return self.time_labels[key]
        try:
            f = float(value)
            if f.is_integer():
                return str(int(f))
        except (TypeError, ValueError):
            pass
        return str(value)

    def add_test(
        self,
        *,
        name: str,
        group_col: str,
        subset: dict[str, Any] | None = None,
        levels: Sequence[Any] | None = None,
        comparisons: Sequence[tuple[Any, Any]] | str | None = "all",
        min_subjects_per_group: int = 2,
        run_kruskal: bool = True,
        run_pairwise: bool = True,
    ) -> "BaselineDistanceAnalysis":
        if group_col not in self.metadata.columns:
            raise KeyError(f"Group column '{group_col}' is not present in metadata")

        subset = deepcopy(subset or {})
        for column in subset:
            if column not in self.metadata.columns:
                raise KeyError(f"Subset column '{column}' is not present in metadata")
        if group_col in subset:
            raise ValueError("Cannot subset on the same column that is being tested")
        if any(t["name"] == name for t in self.tests):
            raise ValueError(f"A test named '{name}' already exists")

        if isinstance(comparisons, str):
            if comparisons != "all":
                raise ValueError("comparisons must be 'all', None, or a sequence of pairs")
            comparison_value = "all"
        elif comparisons is None:
            comparison_value = None
        else:
            comparison_value = [(str(a), str(b)) for a, b in comparisons]

        self.tests.append({
            "name": name,
            "group_col": group_col,
            "subset": subset,
            "levels": [str(x) for x in levels] if levels is not None else None,
            "comparisons": comparison_value,
            "min_subjects_per_group": int(min_subjects_per_group),
            "run_kruskal": bool(run_kruskal),
            "run_pairwise": bool(run_pairwise),
        })
        return self

    def add_tests_from_columns(
        self,
        *,
        test_cols: Iterable[str],
        subset_cols: Iterable[str],
        include_unfiltered: bool = True,
        max_subset_columns: int | None = None,
        min_subjects_per_group: int = 2,
        level_orders: dict[str, Sequence[Any]] | None = None,
        comparisons: dict[str, Sequence[tuple[Any, Any]] | str | None] | None = None,
        skip_duplicates: bool = True,
    ) -> "BaselineDistanceAnalysis":
        for group_col in test_cols:
            subsets = generate_subsets(
                self.metadata,
                subset_cols,
                exclude_cols=[group_col],
                include_unfiltered=include_unfiltered,
                max_subset_columns=max_subset_columns,
            )

            for subset in subsets:
                duplicate = any(
                    t["group_col"] == group_col and t["subset"] == subset
                    for t in self.tests
                )
                if duplicate and skip_duplicates:
                    continue

                desc = subset_description(subset)
                safe = (
                    desc.replace(" & ", "__")
                    .replace("=", "-")
                    .replace(" ", "_")
                    .replace("[", "")
                    .replace("]", "")
                    .replace(",", "-")
                )
                name = f"{group_col}__{safe}"

                self.add_test(
                    name=name,
                    group_col=group_col,
                    subset=subset,
                    levels=(level_orders or {}).get(group_col),
                    comparisons=(comparisons or {}).get(group_col, "all"),
                    min_subjects_per_group=min_subjects_per_group,
                )
        return self

    def build_distance_table(
        self,
        *,
        time_values: Sequence[Any] | None = None,
        include_baseline: bool = True,
        save: bool = True,
        verbose: bool = True,
    ) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        skipped: list[dict[str, Any]] = []

        for metric, dm in self.distances.items():
            md = self.metadata[
                self.metadata[self.sample_col].astype(str).isin(dm.index)
            ].copy()
            md = md.dropna(subset=[self.sample_col, self.subject_col, self.time_col])

            if verbose:
                print(
                    f"{metric}: {len(dm)} distance-matrix samples, "
                    f"{len(md)} metadata rows, {md[self.subject_col].nunique()} subjects"
                )

            for subject, subject_df in md.groupby(self.subject_col, sort=False):
                baseline_rows = subject_df[
                    self._time_mask(subject_df[self.time_col], self.baseline_time)
                ].copy()

                if baseline_rows.empty:
                    message = f"No baseline sample at {self.baseline_time}"
                    skipped.append({"metric": metric, "subject": subject, "reason": message})
                    if self.missing_baseline == "error":
                        raise ValueError(f"{metric}: subject {subject}: {message}")
                    continue

                if len(baseline_rows) > 1:
                    if self.duplicate_baseline == "error":
                        raise ValueError(
                            f"{metric}: subject {subject} has {len(baseline_rows)} baseline samples "
                            f"at time {self.baseline_time}"
                        )
                    baseline_rows = baseline_rows.sort_values(self.sample_col).iloc[[0]]

                baseline_sample = str(baseline_rows.iloc[0][self.sample_col])

                for _, row in subject_df.iterrows():
                    tp = row[self.time_col]
                    is_baseline = bool(self._time_mask(pd.Series([tp]), self.baseline_time).iloc[0])

                    if not include_baseline and is_baseline:
                        continue

                    if time_values is not None:
                        keep = any(
                            bool(self._time_mask(pd.Series([tp]), wanted).iloc[0])
                            for wanted in time_values
                        )
                        if not keep:
                            continue

                    sample_id = str(row[self.sample_col])
                    distance = float(dm.loc[baseline_sample, sample_id])

                    out = row.to_dict()
                    out.update({
                        "metric": metric,
                        "subject": str(subject),
                        "sample_id": sample_id,
                        "baseline_sample": baseline_sample,
                        "timepoint_value": tp,
                        "timepoint_label": self._time_label(tp),
                        "is_baseline": is_baseline,
                        "dist_to_baseline": distance,
                    })
                    rows.append(out)

        table = pd.DataFrame(rows)
        if table.empty:
            raise RuntimeError("No distance-to-baseline rows were generated")

        table["_time_sort"] = pd.to_numeric(table["timepoint_value"], errors="coerce")
        table = table.sort_values(
            ["metric", "subject", "_time_sort", "sample_id"],
            na_position="last",
        ).drop(columns="_time_sort").reset_index(drop=True)

        self.distance_table_ = table
        self.baseline_skipped_ = pd.DataFrame(skipped)

        if save:
            if self.output_dir is None:
                raise ValueError("output_dir must be configured when save=True")
            self.output_dir.mkdir(parents=True, exist_ok=True)
            table.to_csv(
                self.output_dir / "distance_to_baseline_table.tsv",
                sep="\t",
                index=False,
            )
            table.to_csv(
                self.output_dir / "distance_to_baseline_table.csv",
                index=False,
            )
            self.baseline_skipped_.to_csv(
                self.output_dir / "skipped_baselines.tsv",
                sep="\t",
                index=False,
            )

        return table

    def _subject_time_values(self, df: pd.DataFrame, group_col: str) -> pd.DataFrame:
        # One independent value per subject and timepoint. If technical/biological
        # replicate samples exist at a subject-timepoint, mean their distances.
        grouped = (
            df.groupby(
                ["metric", "subject", "timepoint_value", "timepoint_label", group_col],
                dropna=False,
                as_index=False,
            )["dist_to_baseline"]
            .mean()
        )
        return grouped

    def run(
        self,
        *,
        distance_table: pd.DataFrame | None = None,
        include_baseline_tests: bool = False,
        save: bool = True,
        verbose: bool = True,
    ) -> BaselineResults:
        if not self.tests:
            raise ValueError("No tests configured. Use add_test() or add_tests_from_columns().")

        if distance_table is None:
            if hasattr(self, "distance_table_"):
                distance_table = self.distance_table_.copy()
            else:
                distance_table = self.build_distance_table(save=save, verbose=verbose)
        else:
            distance_table = distance_table.copy()

        if not include_baseline_tests:
            test_table = distance_table[~distance_table["is_baseline"].astype(bool)].copy()
        else:
            test_table = distance_table.copy()

        kw_rows: list[dict[str, Any]] = []
        pw_rows: list[dict[str, Any]] = []
        skipped: list[dict[str, Any]] = []

        for metric in list(pd.unique(test_table["metric"])):
            metric_df = test_table[test_table["metric"] == metric].copy()

            for test in self.tests:
                sub = apply_subset(metric_df, test["subset"])
                group_col = test["group_col"]
                sub = sub.dropna(subset=[group_col, "dist_to_baseline", "subject", "timepoint_value"])

                if sub.empty:
                    skipped.append({
                        "metric": metric,
                        "test_name": test["name"],
                        "reason": "no rows after subsetting",
                    })
                    continue

                # A subject should belong to one tested group in this test.
                counts = (
                    sub[["subject", group_col]]
                    .drop_duplicates()
                    .groupby("subject")[group_col]
                    .nunique()
                )
                bad_subjects = counts[counts > 1].index.tolist()
                if bad_subjects:
                    skipped.append({
                        "metric": metric,
                        "test_name": test["name"],
                        "reason": f"group column changes within subject; examples: {bad_subjects[:10]}",
                    })
                    continue

                values_df = self._subject_time_values(sub, group_col)

                if test["levels"] is not None:
                    levels = list(test["levels"])
                else:
                    levels = list(pd.unique(values_df[group_col].dropna().astype(str)))
                    try:
                        levels = sorted(levels)
                    except TypeError:
                        pass

                if test["comparisons"] == "all":
                    comparisons_to_run = list(combinations(levels, 2))
                elif test["comparisons"] is None:
                    comparisons_to_run = []
                else:
                    comparisons_to_run = list(test["comparisons"])

                time_values = list(pd.unique(values_df["timepoint_value"]))
                try:
                    time_values = sorted(time_values, key=lambda x: float(x))
                except (TypeError, ValueError):
                    time_values = sorted(time_values, key=str)

                for tp in time_values:
                    tp_df = values_df[values_df["timepoint_value"].astype(str) == str(tp)].copy()
                    tp_label = str(tp_df["timepoint_label"].iloc[0]) if not tp_df.empty else str(tp)

                    group_arrays: list[np.ndarray] = []
                    group_names: list[str] = []
                    for level in levels:
                        arr = tp_df.loc[
                            tp_df[group_col].astype(str) == str(level),
                            "dist_to_baseline",
                        ].dropna().to_numpy(dtype=float)
                        if len(arr) >= test["min_subjects_per_group"]:
                            group_arrays.append(arr)
                            group_names.append(str(level))

                    base = {
                        "metric": metric,
                        "test_name": test["name"],
                        "group_col": group_col,
                        "subset_description": subset_description(test["subset"]),
                        "timepoint_value": tp,
                        "timepoint_label": tp_label,
                    }

                    if test["run_kruskal"]:
                        if len(group_arrays) >= 2:
                            H, p = kruskal(*group_arrays)
                            kw_rows.append({
                                **base,
                                "groups_tested": " | ".join(group_names),
                                "n_groups": len(group_arrays),
                                "H_statistic": float(H),
                                "p_value": float(p),
                            })
                        else:
                            skipped.append({
                                **base,
                                "reason": "fewer than 2 groups with sufficient subjects for Kruskal-Wallis",
                            })

                    if test["run_pairwise"]:
                        for group_1, group_2 in comparisons_to_run:
                            x = tp_df.loc[
                                tp_df[group_col].astype(str) == str(group_1),
                                "dist_to_baseline",
                            ].dropna().to_numpy(dtype=float)
                            y = tp_df.loc[
                                tp_df[group_col].astype(str) == str(group_2),
                                "dist_to_baseline",
                            ].dropna().to_numpy(dtype=float)

                            if (
                                len(x) < test["min_subjects_per_group"]
                                or len(y) < test["min_subjects_per_group"]
                            ):
                                skipped.append({
                                    **base,
                                    "group_1": group_1,
                                    "group_2": group_2,
                                    "reason": "insufficient subjects for pairwise test",
                                })
                                continue

                            U, p = mannwhitneyu(x, y, alternative="two-sided")
                            pw_rows.append({
                                **base,
                                "group_1": str(group_1),
                                "group_2": str(group_2),
                                "n_1": int(len(x)),
                                "n_2": int(len(y)),
                                **summarize_two_groups(x, y),
                                "U_statistic": float(U),
                                "p_value": float(p),
                            })

        kw = pd.DataFrame(kw_rows)
        pw = pd.DataFrame(pw_rows)
        kw, pw = add_fdr_columns(kw, pw)
        skipped_df = pd.DataFrame(skipped)

        results = BaselineResults(
            distance_table=distance_table,
            kruskal=kw,
            pairwise=pw,
            skipped=skipped_df,
        )
        self.results_ = results

        if save:
            if self.output_dir is None:
                raise ValueError("output_dir must be configured when save=True")
            self.output_dir.mkdir(parents=True, exist_ok=True)
            kw.to_csv(self.output_dir / "kruskal_wallis.tsv", sep="\t", index=False)
            pw.to_csv(self.output_dir / "pairwise_mannwhitney.tsv", sep="\t", index=False)
            skipped_df.to_csv(self.output_dir / "skipped_tests.tsv", sep="\t", index=False)

        if verbose:
            print(f"Kruskal-Wallis rows: {len(kw)}")
            print(f"Pairwise Mann-Whitney rows: {len(pw)}")
            print(f"Skipped entries: {len(skipped_df)}")

        return results

    def plot(
        self,
        *,
        results: BaselineResults | None = None,
        output_dir: str | Path | None = None,
        q_column: str = "q_value_fdr_bh_per_timepoint",
        alpha: float = 0.05,
        show_only_significant: bool = True,
        group_colors: dict[str, str] | None = None,
        group_labels: dict[str, str] | None = None,
        time_order: Sequence[Any] | None = None,
        time_labels: dict[Any, str] | None = None,
        y_lower: float | None = None,
        save_png: bool = True,
        save_pdf: bool = True,
        show: bool = False,
    ):
        from .plotting import plot_all

        if results is None:
            if not hasattr(self, "results_"):
                raise ValueError("Run the analysis first or pass results=...")
            results = self.results_

        if output_dir is None:
            if self.output_dir is None:
                raise ValueError("output_dir must be provided")
            output_dir = self.output_dir / "plots"

        return plot_all(
            analysis=self,
            results=results,
            output_dir=output_dir,
            q_column=q_column,
            alpha=alpha,
            show_only_significant=show_only_significant,
            group_colors=group_colors,
            group_labels=group_labels,
            time_order=time_order,
            time_labels=time_labels,
            y_lower=y_lower,
            save_png=save_png,
            save_pdf=save_pdf,
            show=show,
        )
