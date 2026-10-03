"""Tests of the experiments."""

import numpy as np
import pytest

from wavebreak import analysis
from wavebreak.data import load_transects


@pytest.fixture(scope="module")
def transects():
    return load_transects()


def test_calibrated_amplitude():
    assert analysis.calibrated_amplitude(2.0) == pytest.approx(1.02 * 0.78 * 2.0)


def test_transect_table(transects):
    table = analysis.run_transects(transects)
    assert list(table["transect"]) == ["A", "B", "C", "D"]
    assert table["broke"].all()
    assert table["valid"].all()
    assert analysis.ranking(table) == "C>D>B>A"


def test_no_single_wave_height_fits_all_transects(transects):
    """The four validity windows have no common value: a 'same wave everywhere' test is impossible."""
    windows = [analysis.amplitude_window(t) for t in transects.values()]
    assert max(w[0] for w in windows) > min(w[1] for w in windows)


def test_slope_sweep_is_not_monotonic():
    sweep = analysis.slope_sweep(np.array([4.0, 8.0, 20.0]))
    y = sweep["tube_area_norm"].to_numpy()
    assert y[1] > y[0] and y[1] > y[2]  # rises, then falls


def test_monte_carlo_is_reproducible(transects):
    a = analysis.monte_carlo(transects, 5, seed=1)
    b = analysis.monte_carlo(transects, 5, seed=1)
    c = analysis.monte_carlo(transects, 5, seed=2)
    assert a.equals(b)
    assert not a["reef_slope_deg"].equals(c["reef_slope_deg"])


def test_perturbation_stays_within_rounding_error(transects):
    rng = np.random.default_rng(0)
    t = transects["A"]
    for _ in range(100):
        p = analysis.perturbed(t, rng)
        assert np.all(np.abs(p.depths - t.depths) <= 2 * analysis.DEPTH_ERROR + 1e-12)


def test_keep_valid_samples_drops_whole_sample(transects):
    samples = analysis.monte_carlo(transects, 8, seed=3)
    valid = analysis.keep_valid_samples(samples)
    assert valid["valid"].all()
    assert len(valid) % 4 == 0  # complete samples only (4 transects each)


def test_ranking_frequencies_sum_to_one(transects):
    samples = analysis.keep_valid_samples(analysis.monte_carlo(transects, 10, seed=4))
    freq = analysis.ranking_frequencies(samples)
    assert freq.sum() == pytest.approx(1.0)
    assert len(freq) == 24
