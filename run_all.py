"""Run every experiment and write the tables, figures and key numbers.

    python run_all.py                  # full run (about 2 minutes)
    python run_all.py --quick          # fast check (a few seconds), fewer Monte Carlo samples
    python run_all.py --out my_folder  # choose where the results go (default: results/)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from wavebreak import analysis, plots
from wavebreak.data import load_survey_shares, load_transects
from wavebreak.model import BEACH, FIELD, plane_beach, simulate, survey_profile

SEED = 12345  # fixed: the Monte Carlo gives the same numbers every time
N_SAMPLES = 1000
N_SAMPLES_QUICK = 50


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="fewer simulations")
    parser.add_argument("--out", default="results", help="output folder")
    args = parser.parse_args()

    out = Path(args.out)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)
    key: dict = {"seed": SEED, "quick": args.quick}

    def save_table(frame, name):
        frame.to_csv(out / "tables" / f"{name}.csv", index=False, float_format="%.6g")

    def save_figure(fig, name):
        fig.savefig(out / "figures" / f"{name}.png")
        plots.plt.close(fig)

    # 1. Field data ---------------------------------------------------------------
    transects = load_transects()
    save_figure(plots.plot_field(transects), "field_bathymetry")

    # 2. The original study: one calibrated wave per transect ---------------------
    table = analysis.run_transects(transects)
    save_table(table, "transects")
    save_figure(plots.plot_transects(table), "transects")
    lows, highs = table["window_low"], table["window_high"]
    key["transects"] = {
        "ranking_tube_area_norm": analysis.ranking(table),
        "ranking_jet_reach": analysis.ranking(table, "jet_reach"),
        "all_valid": bool(table["valid"].all()),
        "all_overturned_early": bool(table["overturned_early"].all()),
        "common_amplitude_window_exists": bool(lows.max() < highs.min()),
    }

    # 3. Picture of one breaking wave ---------------------------------------------
    beach = plane_beach(6.0, 10.0)
    recorded = simulate(beach, 2.0, BEACH, record=True)
    save_figure(
        plots.plot_sequence(recorded, beach, "One breaking wave on a 10° plane beach (A = 2 m)"),
        "breaking_sequence",
    )

    # 4. One wave, many slopes (controlled experiment) ----------------------------
    slopes = np.arange(2.0, 31.0, 2.0 if args.quick else 1.0)
    sweep = analysis.slope_sweep(slopes)
    save_table(sweep, "slope_sweep")
    resolution = analysis.resolution_check()
    save_table(resolution, "resolution_check")
    values = resolution["tube_area_norm"].to_numpy(dtype=float)
    band = float((values.max() - values.min()) / 2 / values.mean())
    save_figure(plots.plot_slope_sweep(sweep, band), "slope_sweep")
    best = sweep.loc[sweep["tube_area_norm"].idxmax()]
    key["slope_sweep"] = {
        "best_slope_deg": float(best["slope_deg"]),
        "best_tube_area_norm": float(best["tube_area_norm"]),
        "numerical_band": band,
        "slopes_overturned_early_deg": sweep.loc[sweep["overturned_early"], "slope_deg"].tolist(),
    }

    # 5. Sensitivity to the width of the incident wave ----------------------------
    width = analysis.width_sensitivity(transects)
    save_table(width, "width_sensitivity")
    save_figure(plots.plot_width(width), "width_sensitivity")
    key["width_rankings"] = width.groupby("width")["ranking"].first().to_dict()

    # 6. Monte Carlo on the measurement errors ------------------------------------
    n = N_SAMPLES_QUICK if args.quick else N_SAMPLES
    samples = analysis.monte_carlo(transects, n, SEED)
    valid = analysis.keep_valid_samples(samples)
    freq = analysis.ranking_frequencies(valid)
    save_table(samples, "monte_carlo_samples")
    save_table(freq.reset_index().rename(columns={"index": "ranking"}), "monte_carlo_rankings")
    save_figure(plots.plot_monte_carlo(samples, valid, freq), "monte_carlo")
    slopes_ci = samples.groupby("transect")["reef_slope_deg"].quantile([0.025, 0.975]).unstack()
    last = valid.pivot(index="sample", columns="transect", values="tube_area_norm")
    key["monte_carlo"] = {
        "n_samples": n,
        "n_valid": int(valid["sample"].nunique()),
        "invalid_share_per_transect": (1 - samples.groupby("transect")["valid"].mean()).to_dict(),
        "ranking_share": freq.head(4).to_dict(),
        "share_A_last": float((last.idxmin(axis=1) == "A").mean()),
        "slope_95_interval_deg": {k: [float(a), float(b)] for k, (a, b) in slopes_ci.iterrows()},
    }

    # 7. Survey of surfers --------------------------------------------------------
    shares = load_survey_shares()
    save_figure(plots.plot_survey(shares), "survey")

    (out / "key_numbers.json").write_text(json.dumps(key, indent=2, default=float) + "\n")
    print(f"Done. Results in {out}/ (tables, figures, key_numbers.json).")


if __name__ == "__main__":
    main()
