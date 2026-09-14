import pandas as pd
from baselinebeta.subsets import generate_subsets, apply_subset


def test_generate_subsets_two_binary_columns():
    md = pd.DataFrame({
        "model": ["WT", "WT", "APC", "APC"],
        "sex": ["female", "male", "female", "male"],
    })
    subsets = generate_subsets(md, ["model", "sex"])
    assert {} in subsets
    assert {"model": "WT"} in subsets
    assert {"sex": "female"} in subsets
    assert {"model": "WT", "sex": "female"} in subsets


def test_apply_subset():
    df = pd.DataFrame({"model": ["WT", "APC"], "sex": ["female", "female"]})
    out = apply_subset(df, {"model": "WT", "sex": "female"})
    assert len(out) == 1
