"""Tests of the data loading and checks."""

import numpy as np
import pandas as pd
import pytest

from wavebreak.data import (
    check_bathymetry,
    load_survey_shares,
    load_transects,
    reef_face_slope_deg,
)


def test_four_transects_of_four_points():
    t = load_transects()
    assert list(t) == ["A", "B", "C", "D"]
    assert all(len(x.depths) == 4 for x in t.values())
    assert all(x.spacing == 12.0 for x in t.values())


def test_depths_are_measured_minus_tide():
    table = pd.read_csv("data/bathymetry.csv")
    recomputed = (table["measured_depth_m"] - table["tide_correction_m"]).round(1)
    assert np.allclose(recomputed, table["corrected_depth_m"])


def test_reef_face_slopes():
    t = load_transects()
    slopes = {k: reef_face_slope_deg(v.depths, v.spacing) for k, v in t.items()}
    assert slopes["A"] == pytest.approx(0.95, abs=0.01)
    assert slopes["B"] == pytest.approx(5.24, abs=0.01)
    assert slopes["C"] == pytest.approx(5.24, abs=0.01)  # B and C have exactly the same slope
    assert slopes["D"] == pytest.approx(0.48, abs=0.01)


def test_reef_depth_is_the_shallowest_point():
    t = load_transects()
    assert t["A"].reef_depth == pytest.approx(1.5)
    assert t["B"].reef_depth == pytest.approx(2.5)


def test_check_detects_inconsistent_table():
    table = pd.read_csv("data/bathymetry.csv")
    check_bathymetry(table)  # the real table is fine
    table.loc[0, "corrected_depth_m"] += 0.5
    with pytest.raises(ValueError, match="minus tide"):
        check_bathymetry(table)


def test_survey_shares():
    s = load_survey_shares()
    sw = s[(s["question_id"] == "preferred_swell_direction") & (s["option"] == "SW")]
    assert int(sw["count"].iloc[0]) == 26
    assert float(sw["share"].iloc[0]) == pytest.approx(26 / 30)
