from baselinebeta import BaselineDistanceAnalysis

BASE = "/vol/jlab/MicrobiomeAnalyses/Projects/JansenVanVuuren_Coldshield/Karl"

analysis = BaselineDistanceAnalysis(
    metadata=f"{BASE}/metadata_coldshield_behaviour.txt",
    distances={
        "unweighted_unifrac": f"{BASE}/beta_diversity/unweighted_unifrac/unweighted-unifrac-distance-matrix.qza",
    },
    sample_col="sample_name",
    subject_col="host_subject_id",
    time_col="timepoint_numeric",
    baseline_time=-1,
    output_dir=f"{BASE}/baselinebeta_results",
    metric_titles={"unweighted_unifrac": "Unweighted UniFrac"},
    time_labels={-1: "Baseline", 1: "Day 1", 3: "Day 3", 7: "Day 7"},
    missing_baseline="skip",
    duplicate_baseline="error",
)

analysis.add_test(
    name="treatment",
    group_col="description_of_treatment",
    levels=["cold", "cold_radiation", "warm", "warm_radiation"],
    comparisons=[
        ("cold", "cold_radiation"),
        ("cold", "warm"),
        ("cold_radiation", "warm_radiation"),
        ("warm", "warm_radiation"),
    ],
)

distance_table = analysis.build_distance_table(time_values=[-1, 1, 3, 7])
results = analysis.run(distance_table=distance_table)

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
    show_only_significant=True,
)
