import pandas as pd
from baselinebeta import BaselineDistanceAnalysis


def test_distance_to_baseline_table():
    ids = ["A0", "A1", "B0", "B1"]
    dm = pd.DataFrame(
        [
            [0.0, 0.2, 0.4, 0.5],
            [0.2, 0.0, 0.3, 0.4],
            [0.4, 0.3, 0.0, 0.6],
            [0.5, 0.4, 0.6, 0.0],
        ],
        index=ids,
        columns=ids,
    )
    md = pd.DataFrame({
        "sample": ids,
        "subject": ["A", "A", "B", "B"],
        "time": [0, 1, 0, 1],
        "group": ["x", "x", "y", "y"],
    })

    analysis = BaselineDistanceAnalysis(
        metadata=md,
        distances={"d": dm},
        sample_col="sample",
        subject_col="subject",
        time_col="time",
        baseline_time=0,
        duplicate_baseline="error",
    )
    table = analysis.build_distance_table(save=False, verbose=False)
    a1 = table.loc[table["sample_id"] == "A1", "dist_to_baseline"].iloc[0]
    b1 = table.loc[table["sample_id"] == "B1", "dist_to_baseline"].iloc[0]
    assert a1 == 0.2
    assert b1 == 0.6
