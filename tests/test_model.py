"""Tests of the physical model: do the simulated numbers match what we can compute by hand?"""

import numpy as np
import pytest

from wavebreak.model import (
    BEACH,
    FIELD,
    G,
    Params,
    Result,
    plane_beach,
    simulate,
    starts_breaking,
    survey_profile,
)


def flat(h):
    return lambda x: h + 0.0 * np.asarray(x)


# ------------------------------------------------------------------ the sea bed


def test_plane_beach_depth():
    depth = plane_beach(6.0, 10.0)
    assert depth(np.array([-5.0, 0.0]))[0] == pytest.approx(6.0)  # flat before the toe
    assert depth(np.array([10.0]))[0] == pytest.approx(6.0 - 10.0 * np.tan(np.radians(10.0)))


def test_survey_profile_interpolates_and_stays_flat_outside():
    depth = survey_profile(np.array([2.0, 3.0, 1.0]), spacing=10.0)
    assert depth(np.array([5.0]))[0] == pytest.approx(2.5)
    assert depth(np.array([-20.0]))[0] == pytest.approx(2.0)
    assert depth(np.array([99.0]))[0] == pytest.approx(1.0)


# ------------------------------------------------------------ simple physics


def test_small_wave_travels_at_sqrt_g_h():
    """A tiny wave on a flat bed moves at the shallow-water speed sqrt(g h)."""
    h, steps = 2.0, 300
    p = Params(nonlinear_from=None, max_steps=steps)
    res = simulate(flat(h), 1e-4, p, record=True)
    x_crest = [x[np.argmax(e)] for x, e in zip(res.x_history, res.eta_history)]
    speed = (x_crest[steps] - x_crest[0]) / (steps * p.dt)
    assert speed == pytest.approx(np.sqrt(G * h), rel=0.01)


def test_flat_bed_breaking_time_matches_theory():
    """On a flat bed, the front overturns at a time known from the method of characteristics.

    Each point moves at c(eta0) = sqrt(g (h + eta0)). The first time two points meet is
        t = 1 / max over x0 of ( - d c / d x0 ),   with  d c / d x0 = g / (2 c) * d eta0 / d x0.
    """
    h, amplitude, sigma = 2.0, 0.5, 3.0
    p = Params(width=sigma, nonlinear_from=None, lip_threshold=1e-3, n_markers=1800, x_min=-60.0)
    x0 = np.linspace(-60.0, 20.0, 400001)
    eta0 = amplitude * np.exp(-((x0 + p.crest_offset) ** 2) / (2 * sigma**2))
    c0 = np.sqrt(G * (h + eta0))
    dc_dx = np.gradient(c0, x0)
    t_theory = 1.0 / np.max(-dc_dx)

    res = simulate(flat(h), amplitude, p)
    assert res.first_crossing_step is not None
    t_model = res.time(res.first_crossing_step)
    assert t_model == pytest.approx(t_theory, rel=0.005)  # observed: 0.04 %


def test_lip_flight_is_a_projectile():
    """The distance flown by the closing lip is u * sqrt(2 H / g) (free fall from height H)."""
    depth = plane_beach(6.0, 10.0)
    res = simulate(depth, 2.0, BEACH)
    assert res.broke
    h_local = float(depth(res.release_x))
    speed = (1.0 + BEACH.jet_boost) * np.sqrt(G * (h_local + res.release_height))
    expected = speed * np.sqrt(2.0 * res.release_height / G)
    assert res.lip_throw == pytest.approx(expected, rel=0.01)  # observed: 0.4 %


def test_tube_area_formula():
    """Area under a parabola: 2/3 * height * length."""
    r = Result(amplitude=2.0, dt=0.01, release_height=1.5, release_x=10.0, closure_x=13.0,
               closure_step=100, onset_x=8.0)
    assert r.lip_throw == pytest.approx(3.0)
    assert r.tube_area == pytest.approx(2.0 / 3.0 * 1.5 * 3.0)
    assert r.tube_area_normalized == pytest.approx(r.tube_area / 4.0)
    assert r.jet_reach == pytest.approx(5.0)


# --------------------------------------------------------------- breaking logic


def test_small_wave_over_deep_flat_bed_does_not_break():
    """eta / h = 0.3 / 6 is far below gamma = 0.78 and the depth never changes."""
    res = simulate(flat(6.0), 0.3, BEACH)
    assert not res.broke
    assert res.onset_step is None


def test_big_wave_breaks_and_onset_is_before_closure():
    res = simulate(plane_beach(6.0, 10.0), 2.0, BEACH)
    assert res.broke
    assert res.onset_step < res.closure_step
    assert res.closure_x > res.onset_x
    assert res.tube_area > 0


def test_starts_breaking_flag():
    depth = flat(2.0)
    assert starts_breaking(depth, 0.78 * 2.0 + 0.01, FIELD)
    assert not starts_breaking(depth, 0.78 * 2.0 - 0.01, FIELD)


def test_onset_respects_gamma_criterion():
    """At onset the highest point satisfies eta >= gamma * local depth."""
    depth = plane_beach(6.0, 10.0)
    res = simulate(depth, 2.0, BEACH, record=True)
    k = res.onset_step + 1  # history index of the state computed at onset_step
    top = int(np.argmax(res.eta_history[k]))
    assert res.eta_history[k][top] >= BEACH.gamma * float(depth(res.x_history[k][top]))


# ----------------------------------------------------- regression on known values


def test_matches_the_original_study_values():
    """Values produced by the first (more complex) version of this project and by the
    original TIPE script logic; the simplification must not have changed them."""
    depth = survey_profile(np.array([2.8, 3.2, 2.1, 2.3]), 12.0)  # transect C
    res = simulate(depth, 1.02 * 0.78 * 2.1)
    assert res.onset_x == pytest.approx(23.5522, abs=1e-3)
    assert res.closure_x == pytest.approx(27.8986, abs=1e-3)
    assert res.jet_reach == pytest.approx(4.3464, abs=1e-3)
    assert res.tube_area == pytest.approx(1.2998, abs=1e-3)

    res7 = simulate(plane_beach(6.0, 7.0), 2.0, BEACH)
    assert res7.tube_area_normalized == pytest.approx(0.5642, abs=1e-3)


def test_original_results_table_is_reproduced():
    """Results table of the original TIPE script (2 decimals): onset position and jet reach."""
    original = {  # transect: (depths, onset_x, jet_reach)
        "A": ([2.1, 1.9, 1.7, 1.5], 34.2, 0.50),
        "B": ([2.6, 3.0, 3.6, 2.5], 35.5, 2.05),
        "C": ([2.8, 3.2, 2.1, 2.3], 23.6, 4.35),
        "D": ([2.4, 2.0, 1.9, 2.1], 19.5, 4.08),
    }
    for name, (depths, onset_x, jet_reach) in original.items():
        depths = np.array(depths)
        res = simulate(survey_profile(depths, 12.0), 1.02 * 0.78 * depths.min())
        assert res.onset_x == pytest.approx(onset_x, abs=0.05), name
        assert res.jet_reach == pytest.approx(jet_reach, abs=0.005), name
