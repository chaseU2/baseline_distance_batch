from __future__ import annotations

import argparse
from pathlib import Path

from .analysis import BaselineDistanceAnalysis
from .config import load_yaml


def build_from_config(cfg: dict) -> BaselineDistanceAnalysis:
    analysis = BaselineDistanceAnalysis(
        metadata=cfg["metadata"],
        distances=cfg["distances"],
        sample_col=cfg["sample_col"],
        subject_col=cfg["subject_col"],
        time_col=cfg["time_col"],
        baseline_time=cfg["baseline_time"],
        output_dir=cfg.get("output_dir"),
        metric_titles=cfg.get("metric_titles"),
        time_labels=cfg.get("time_labels"),
        duplicate_baseline=cfg.get("duplicate_baseline", "error"),
        missing_baseline=cfg.get("missing_baseline", "skip"),
    )

    for test in cfg.get("tests", []):
        analysis.add_test(
            name=test["name"],
            group_col=test["group_col"],
            subset=test.get("subset"),
            levels=test.get("levels"),
            comparisons=test.get("comparisons", "all"),
            min_subjects_per_group=test.get("min_subjects_per_group", 2),
            run_kruskal=test.get("run_kruskal", True),
            run_pairwise=test.get("run_pairwise", True),
        )

    auto = cfg.get("auto_tests")
    if auto:
        analysis.add_tests_from_columns(
            test_cols=auto["test_cols"],
            subset_cols=auto.get("subset_cols", []),
            include_unfiltered=auto.get("include_unfiltered", True),
            max_subset_columns=auto.get("max_subset_columns"),
            min_subjects_per_group=auto.get("min_subjects_per_group", 2),
            level_orders=auto.get("level_orders"),
            comparisons=auto.get("comparisons"),
        )

    return analysis


def main() -> None:
    parser = argparse.ArgumentParser(description="Distance-to-baseline beta-diversity analysis")
    parser.add_argument("config", help="YAML configuration file")
    parser.add_argument("--no-plots", action="store_true", help="Do not generate plots")
    args = parser.parse_args()

    cfg = load_yaml(args.config)
    analysis = build_from_config(cfg)
    table = analysis.build_distance_table(
        time_values=cfg.get("time_values"),
        include_baseline=True,
        save=True,
    )
    results = analysis.run(distance_table=table, save=True)

    if not args.no_plots:
        plot_cfg = cfg.get("plot", {})
        analysis.plot(
            results=results,
            q_column=plot_cfg.get("q_column", "q_value_fdr_bh_per_timepoint"),
            alpha=plot_cfg.get("alpha", 0.05),
            show_only_significant=plot_cfg.get("show_only_significant", True),
            group_colors=plot_cfg.get("group_colors"),
            group_labels=plot_cfg.get("group_labels"),
            time_order=plot_cfg.get("time_order"),
            time_labels=plot_cfg.get("time_labels"),
            y_lower=plot_cfg.get("y_lower"),
            show=False,
        )


if __name__ == "__main__":
    main()
