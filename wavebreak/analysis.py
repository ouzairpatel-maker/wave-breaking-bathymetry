"""The experiments: each function runs the model several times and returns a table.

All tables are pandas DataFrames; ``run_all.py`` saves them as CSV files.
"""

from __future__ import annotations

from dataclasses import replace
from itertools import permutations

import numpy as np
import pandas as pd

from wavebreak.data import Transect, reef_face_slope_deg
from wavebreak.model import (
    BEACH,
    FIELD,
    Params,
    Result,
    plane_beach,
    simulate,
    starts_breaking,
    survey_profile,
)

# The original study chose each wave so that it just reaches the breaking criterion
# above the reef: A = MARGIN * gamma * (reef depth).
AMPLITUDE_MARGIN = 1.02


def calibrated_amplitude(reef_depth: float, params: Params = FIELD) -> float:
    """Wave height used in the original study for a reef of the given depth [m]."""
    return AMPLITUDE_MARGIN * params.gamma * reef_depth


def amplitude_window(transect: Transect, params: Params = FIELD) -> tuple[float, float]:
    """Wave heights for which the model is usable on this transect.

    Below the lower bound the wave never reaches eta/h = gamma. At or above the
    upper bound the wave already satisfies the criterion offshore, before the run.
    """
    return params.gamma * transect.reef_depth, params.gamma * transect.offshore_depth


def summarize(result: Result) -> dict:
    """The numbers we keep from one run."""
    return {
        "broke": result.broke,
        "overturned_early": result.overturned_early,
        "onset_x": result.onset_x,
        "jet_reach": result.jet_reach,
        "lip_throw": result.lip_throw,
        "tube_area": result.tube_area,
        "tube_area_norm": result.tube_area_normalized,
    }


# ------------------------------------------------------------------ field transects


def run_transects(transects: dict[str, Transect], params: Params = FIELD) -> pd.DataFrame:
    """Original study: each transect gets its own calibrated wave."""
    rows = []
    for name, t in transects.items():
        depth = survey_profile(t.depths, t.spacing)
        amplitude = calibrated_amplitude(t.reef_depth, params)
        result = simulate(depth, amplitude, params)
        low, high = amplitude_window(t, params)
        rows.append(
            {
                "transect": name,
                "swell": t.swell_direction,
                "offshore_depth": t.offshore_depth,
                "reef_depth": t.reef_depth,
                "reef_slope_deg": reef_face_slope_deg(t.depths, t.spacing),
                "amplitude": amplitude,
                "window_low": low,
                "window_high": high,
                "valid": not starts_breaking(depth, amplitude, params),
                **summarize(result),
            }
        )
    return pd.DataFrame(rows)


def ranking(table: pd.DataFrame, metric: str = "tube_area_norm") -> str:
    """Transect names from the largest to the smallest metric, e.g. 'C>D>B>A'.

    A wave that did not break counts as the smallest.
    """
    values = table.set_index("transect")[metric].astype(float).fillna(-np.inf)
    return ">".join(values.sort_values(ascending=False, kind="stable").index)


# ---------------------------------------------------------------- one wave, many slopes


def slope_sweep(
    slopes_deg: np.ndarray,
    offshore_depth: float = 6.0,
    amplitude: float = 2.0,
    params: Params = BEACH,
) -> pd.DataFrame:
    """Same wave on plane beaches of different slopes: the controlled experiment."""
    rows = []
    for slope in slopes_deg:
        result = simulate(plane_beach(offshore_depth, float(slope)), amplitude, params)
        rows.append({"slope_deg": float(slope), **summarize(result)})
    return pd.DataFrame(rows)


def resolution_check(slope_deg: float = 10.0) -> pd.DataFrame:
    """Run one case with finer and finer grids: how much does the answer move?"""
    rows = []
    for factor in (1, 2, 4):
        params = replace(
            BEACH,
            n_markers=BEACH.n_markers * factor,
            dt=BEACH.dt / factor,
            max_steps=BEACH.max_steps * factor,
        )
        result = simulate(plane_beach(6.0, slope_deg), 2.0, params)
        rows.append(
            {"n_markers": params.n_markers, "dt": params.dt, **summarize(result)}
        )
    return pd.DataFrame(rows)


# ------------------------------------------------------------ sensitivity to wave width


def width_sensitivity(
    transects: dict[str, Transect], widths: tuple[float, ...] = (2.0, 3.0, 4.0, 5.0, 6.0)
) -> pd.DataFrame:
    """Repeat the transect study with other widths of the incident wave."""
    rows = []
    for width in widths:
        table = run_transects(transects, replace(FIELD, width=width))
        order = ranking(table)
        for _, r in table.iterrows():
            rows.append(
                {
                    "width": width,
                    "transect": r["transect"],
                    "tube_area_norm": r["tube_area_norm"],
                    "valid": r["valid"],
                    "ranking": order,
                }
            )
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------ Monte Carlo

DEPTH_ERROR = 0.05  # depths and tide corrections are rounded to 0.1 m -> +/- 0.05 m


def perturbed(t: Transect, rng: np.random.Generator) -> Transect:
    """One random survey compatible with the recorded one.

    Each depth gets its own rounding error; the tide correction error is the same
    for the whole transect (it was read once).
    """
    point_errors = rng.uniform(-DEPTH_ERROR, DEPTH_ERROR, t.depths.size)
    tide_error = rng.uniform(-DEPTH_ERROR, DEPTH_ERROR)
    return replace(t, depths=t.depths + point_errors + tide_error)


def monte_carlo(
    transects: dict[str, Transect], n_samples: int, seed: int, params: Params = FIELD
) -> pd.DataFrame:
    """Repeat the transect study on ``n_samples`` randomly perturbed surveys."""
    rng = np.random.default_rng(seed)
    rows = []
    for sample in range(n_samples):
        survey = {name: perturbed(t, rng) for name, t in transects.items()}
        for name, t in survey.items():
            depth = survey_profile(t.depths, t.spacing)
            amplitude = calibrated_amplitude(t.reef_depth, params)
            result = simulate(depth, amplitude, params)
            rows.append(
                {
                    "sample": sample,
                    "transect": name,
                    "reef_slope_deg": reef_face_slope_deg(t.depths, t.spacing),
                    "amplitude": amplitude,
                    "valid": not starts_breaking(depth, amplitude, params),
                    **summarize(result),
                }
            )
    return pd.DataFrame(rows)


def keep_valid_samples(samples: pd.DataFrame) -> pd.DataFrame:
    """Keep the samples where the wave is usable on all four transects."""
    all_valid = samples.groupby("sample")["valid"].transform("all")
    return samples[all_valid]


def ranking_frequencies(samples: pd.DataFrame, metric: str = "tube_area_norm") -> pd.Series:
    """How often each full ranking (e.g. 'C>D>B>A') occurs, most frequent first."""
    wide = samples.pivot(index="sample", columns="transect", values=metric).astype(float)
    wide = wide.fillna(-np.inf)
    names = list(wide.columns)
    orders = np.argsort(-wide.to_numpy(), axis=1, kind="stable")
    labels = pd.Series([">".join(names[i] for i in row) for row in orders])
    every = [">".join(p) for p in permutations(names)]
    freq = labels.value_counts(normalize=True).reindex(every, fill_value=0.0)
    return freq.sort_values(ascending=False).rename("frequency")
