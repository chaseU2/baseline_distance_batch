from __future__ import annotations

import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def _validate_distance_matrix(df: pd.DataFrame, source: str = "distance matrix") -> pd.DataFrame:
    df = df.copy()
    df.index = df.index.astype(str).str.strip()
    df.columns = df.columns.astype(str).str.strip()
    df = df.astype(float)

    if df.shape[0] != df.shape[1]:
        raise ValueError(f"{source}: distance matrix is not square")

    if set(df.index) != set(df.columns):
        raise ValueError(f"{source}: row and column IDs differ")

    if list(df.index) != list(df.columns):
        df = df.loc[df.index, df.index]

    if not np.allclose(df.values, df.values.T, atol=1e-10, equal_nan=True):
        raise ValueError(f"{source}: distance matrix is not symmetric")

    if not np.allclose(np.diag(df.values), 0.0, atol=1e-10, equal_nan=False):
        raise ValueError(f"{source}: diagonal is not zero")

    return df


def _load_qza(path: Path) -> pd.DataFrame:
    with zipfile.ZipFile(path, "r") as zf:
        candidates = [
            name for name in zf.namelist()
            if name.endswith("/data/distance-matrix.tsv")
        ]
        if not candidates:
            candidates = [
                name for name in zf.namelist()
                if "/data/" in name and name.endswith(".tsv")
            ]
        if len(candidates) != 1:
            raise FileNotFoundError(
                f"Expected exactly one distance-matrix TSV inside {path}; found {len(candidates)}"
            )
        raw = zf.read(candidates[0])

    df = pd.read_csv(io.BytesIO(raw), sep="\t", index_col=0)
    return _validate_distance_matrix(df, str(path))


def load_distance_matrix(source: str | Path | pd.DataFrame) -> pd.DataFrame:
    """Load a square distance matrix from QIIME 2 QZA, TSV/CSV, or DataFrame."""
    if isinstance(source, pd.DataFrame):
        return _validate_distance_matrix(source, "DataFrame")

    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(path)

    suffix = path.suffix.lower()
    if suffix == ".qza":
        return _load_qza(path)
    if suffix in {".tsv", ".txt"}:
        df = pd.read_csv(path, sep="\t", index_col=0)
        return _validate_distance_matrix(df, str(path))
    if suffix == ".csv":
        df = pd.read_csv(path, index_col=0)
        return _validate_distance_matrix(df, str(path))

    raise ValueError(
        f"Unsupported distance-matrix format '{suffix}'. Use .qza, .tsv, .txt, .csv, or a DataFrame."
    )
