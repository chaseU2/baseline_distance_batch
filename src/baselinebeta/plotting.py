from __future__ import annotations

from itertools import combinations
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from .subsets import apply_subset


def significance_stars(q_value: float) -> str:
    if pd.isna(q_value):
        return ""
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return "n.s."


def add_significance_bracket(ax, x1, x2, y, height, label):
    ax.plot(
        [x1, x1, x2, x2],
        [y, y + height, y + height, y],
        color="black",
        linewidth=1.1,
        clip_on=False,
    )
    ax.text(
        (x1 + x2) / 2,
        y + height,
        label,
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )


def _sort_time_values(values: list[Any]) -> list[Any]:
    try:
        return sorted(values, key=lambda x: float(x))
    except (TypeError, ValueError):
        return sorted(values, key=str)


def plot_all(
    *,
    analysis,
    results,
    output_dir: str | Path,
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
) -> list[Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    distance_table = results.distance_table.copy()
    pairwise = results.pairwise.copy()

    if q_column not in pairwise.columns and not pairwise.empty:
        raise KeyError(f"'{q_column}' not found in pairwise result table")

    sns.set_theme(context="notebook", style="whitegrid")
    written: list[Path] = []

    for metric in pd.unique(distance_table["metric"]):
        metric_df = distance_table[
            (distance_table["metric"] == metric)
            & (~distance_table["is_baseline"].astype(bool))
        ].copy()

        for test in analysis.tests:
            group_col = test["group_col"]
            plot_df = apply_subset(metric_df, test["subset"])
            plot_df = plot_df.dropna(subset=[group_col, "dist_to_baseline", "subject"])
            if plot_df.empty:
                continue

            # Aggregate replicate samples to one subject value per timepoint.
            plot_df = (
                plot_df.groupby(
                    ["subject", "timepoint_value", "timepoint_label", group_col],
                    dropna=False,
                    as_index=False,
                )["dist_to_baseline"]
                .mean()
            )

            if test["levels"] is not None:
                levels = [str(x) for x in test["levels"]]
            else:
                levels = list(pd.unique(plot_df[group_col].dropna().astype(str)))
                try:
                    levels = sorted(levels)
                except TypeError:
                    pass

            if len(levels) < 2:
                continue

            plot_df = plot_df[plot_df[group_col].astype(str).isin(levels)].copy()
            plot_df[group_col] = pd.Categorical(
                plot_df[group_col].astype(str),
                categories=levels,
                ordered=True,
            )

            if time_order is None:
                times = _sort_time_values(list(pd.unique(plot_df["timepoint_value"])))
            else:
                times = list(time_order)

            # Keep only times with data.
            times = [
                tp for tp in times
                if (plot_df["timepoint_value"].astype(str) == str(tp)).any()
            ]
            if not times:
                continue

            labels_map = {str(k): str(v) for k, v in (time_labels or {}).items()}
            time_titles = []
            for tp in times:
                if str(tp) in labels_map:
                    time_titles.append(labels_map[str(tp)])
                else:
                    rows = plot_df[plot_df["timepoint_value"].astype(str) == str(tp)]
                    if not rows.empty:
                        time_titles.append(str(rows["timepoint_label"].iloc[0]))
                    else:
                        time_titles.append(str(tp))

            palette = None
            if group_colors is not None:
                palette = {
                    level: group_colors.get(level)
                    for level in levels
                    if level in group_colors
                }
                if len(palette) != len(levels):
                    palette = None

            display_labels = {
                level: (group_labels or {}).get(level, level.replace("_", " ").title())
                for level in levels
            }

            sig_df = pairwise[
                (pairwise["metric"] == metric)
                & (pairwise["test_name"] == test["name"])
            ].copy()
            if not sig_df.empty:
                sig_df[q_column] = pd.to_numeric(sig_df[q_column], errors="coerce")

            if test["comparisons"] == "all":
                comparison_order = list(combinations(levels, 2))
            elif test["comparisons"] is None:
                comparison_order = []
            else:
                comparison_order = list(test["comparisons"])
            pair_rank = {tuple(map(str, pair)): i for i, pair in enumerate(comparison_order)}
            if not sig_df.empty:
                sig_df["_pair_order"] = [
                    pair_rank.get((str(a), str(b)), 999)
                    for a, b in zip(sig_df["group_1"], sig_df["group_2"])
                ]

            y_values = plot_df["dist_to_baseline"].dropna()
            if y_values.empty:
                continue
            data_min = float(y_values.min())
            data_max = float(y_values.max())
            data_range = data_max - data_min
            if not np.isfinite(data_range) or data_range == 0:
                data_range = max(abs(data_max), 1.0)

            bracket_step = 0.10 * data_range
            bracket_height = 0.025 * data_range
            bracket_start = data_max + 0.07 * data_range

            max_brackets = 0
            for tp in times:
                tp_tests = sig_df[sig_df["timepoint_value"].astype(str) == str(tp)].copy()
                if show_only_significant and not tp_tests.empty:
                    tp_tests = tp_tests[tp_tests[q_column] < alpha]
                max_brackets = max(max_brackets, len(tp_tests))

            upper = bracket_start + max(max_brackets, 1) * bracket_step + 0.10 * data_range
            lower = y_lower if y_lower is not None else max(0.0, data_min - 0.05 * data_range)

            fig, axes = plt.subplots(
                1,
                len(times),
                figsize=(max(5.0, 4.5 * len(times)), 5.8),
                sharey=True,
                squeeze=False,
            )
            axes = axes[0]

            for ax, tp, tp_title in zip(axes, times, time_titles):
                sub = plot_df[plot_df["timepoint_value"].astype(str) == str(tp)].copy()

                sns.boxplot(
                    data=sub,
                    x=group_col,
                    y="dist_to_baseline",
                    order=levels,
                    hue=group_col,
                    hue_order=levels,
                    palette=palette,
                    dodge=False,
                    legend=False,
                    width=0.62,
                    showfliers=False,
                    linewidth=1.3,
                    medianprops={"color": "black", "linewidth": 1.6},
                    whiskerprops={"color": "black", "linewidth": 1.1},
                    capprops={"color": "black", "linewidth": 1.1},
                    boxprops={"edgecolor": "black"},
                    ax=ax,
                )
                for patch in ax.patches:
                    patch.set_alpha(0.65)

                sns.stripplot(
                    data=sub,
                    x=group_col,
                    y="dist_to_baseline",
                    order=levels,
                    hue=group_col,
                    hue_order=levels,
                    palette=palette,
                    dodge=False,
                    legend=False,
                    jitter=0.18,
                    size=6,
                    alpha=0.90,
                    edgecolor="black",
                    linewidth=0.7,
                    zorder=3,
                    ax=ax,
                )

                counts = (
                    sub.groupby(group_col, observed=False)
                    .size()
                    .reindex(levels, fill_value=0)
                )
                tick_labels = [
                    f"{display_labels[level]}\n$n$ = {int(counts[level])}"
                    for level in levels
                ]
                ax.set_xticks(range(len(levels)))
                ax.set_xticklabels(tick_labels, fontsize=9)

                tp_tests = sig_df[sig_df["timepoint_value"].astype(str) == str(tp)].copy()
                if show_only_significant and not tp_tests.empty:
                    tp_tests = tp_tests[tp_tests[q_column] < alpha].copy()
                if not tp_tests.empty:
                    tp_tests = tp_tests.sort_values("_pair_order")

                for bracket_number, (_, row) in enumerate(tp_tests.iterrows()):
                    g1 = str(row["group_1"])
                    g2 = str(row["group_2"])
                    if g1 not in levels or g2 not in levels:
                        continue
                    x1 = levels.index(g1)
                    x2 = levels.index(g2)
                    y = bracket_start + bracket_number * bracket_step
                    add_significance_bracket(
                        ax,
                        x1,
                        x2,
                        y,
                        bracket_height,
                        significance_stars(row[q_column]),
                    )

                if show_only_significant and len(tp_tests) == 0:
                    ax.text(
                        0.5,
                        0.97,
                        "No FDR-significant\npairwise differences",
                        transform=ax.transAxes,
                        ha="center",
                        va="top",
                        fontsize=9,
                        color="dimgray",
                    )

                ax.set_title(tp_title, fontsize=13, fontweight="bold")
                ax.set_xlabel("")
                ax.set_ylim(lower, upper)
                ax.grid(axis="x", visible=False)
                ax.grid(axis="y", linestyle="--", linewidth=0.6, alpha=0.5)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)

            axes[0].set_ylabel("Distance to individual baseline", fontsize=11)
            for ax in axes[1:]:
                ax.set_ylabel("")

            metric_title = analysis.metric_titles.get(
                metric,
                metric.replace("_", " ").title(),
            )
            fig.suptitle(
                f"{metric_title} — {test['name']}",
                fontsize=15,
                fontweight="bold",
                y=1.02,
            )
            fig.text(
                0.5,
                -0.02,
                (
                    "Points represent individual subjects; boxes show the interquartile range and median. "
                    "Significance is based on two-sided Mann–Whitney U tests."
                ),
                ha="center",
                fontsize=9,
            )
            plt.tight_layout()

            stem = f"{metric}__{test['name']}__distance_to_baseline"
            if save_png:
                path = output_dir / f"{stem}.png"
                fig.savefig(path, dpi=300, bbox_inches="tight")
                written.append(path)
            if save_pdf:
                path = output_dir / f"{stem}.pdf"
                fig.savefig(path, bbox_inches="tight")
                written.append(path)

            if show:
                plt.show()
            else:
                plt.close(fig)

    return written
