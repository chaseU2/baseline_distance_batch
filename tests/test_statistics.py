import numpy as np
from baselinebeta.statistics import cliffs_delta


def test_cliffs_delta_direction():
    x = np.array([1, 2, 3])
    y = np.array([4, 5, 6])
    assert cliffs_delta(x, y) == 1.0
    assert cliffs_delta(y, x) == -1.0
