"""All the figures. Each function returns a matplotlib Figure."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # draw to files, no window needed

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from wavebreak.data import Transect
from wavebreak.model import DepthFunction, Result

# One colour per transect, readable by colour-blind people (Okabe-Ito palette).
COLORS = {"A": "#0072B2", "B": "#D55E00", "C": "#009E73", "D": "#CC79A7"}

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "figure.dpi": 100,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    }
)


def plot_sequence(result: Result, depth: DepthFunction, title: str) -> Figure:
    """Six snapshots of one run: the wave steepens, the lip is thrown, the tube closes."""
    assert result.onset_step is not None and result.closure_step is not None
    # x_history[0] is the initial state; x_history[k] is the state computed at step k - 1.
    onset, closure = result.onset_step + 1, result.closure_step + 1
    frames = [
        0,
        onset // 2,
        onset,
        onset + (closure - onset) // 3,
        onset + 2 * (closure - onset) // 3,
        closure,
    ]
    x_bed = np.linspace(-30, 45, 400)
    bed = -depth(x_bed)

    fig, axes = plt.subplots(2, 3, figsize=(12, 5.5), sharex=True, sharey=True)
    for ax, k in zip(axes.ravel(), frames):
        ax.fill_between(x_bed, bed, -8, color="#c9b27c")  # sea bed
        ax.fill_between(x_bed, bed, 0, color="#cfe5f3")  # still water
        ax.plot(result.x_history[k], result.eta_history[k], color="#0b3d91", lw=1.6)
        if k == closure:
            ax.plot(result.closure_x, 0, "o", color="#c0392b", ms=6, label="tube closes")
            ax.legend(loc="upper left", frameon=False)
        label = "start" if k == 0 else f"t = {k * result.dt:.2f} s"
        if k == onset:
            label += "  (breaking onset)"
        ax.set_title(label, fontsize=10)
        ax.set_xlim(-30, 40)
        ax.set_ylim(-7, 4)
        ax.grid(False)
    for ax in axes[1]:
        ax.set_xlabel("x [m]")
    for ax in axes[:, 0]:
        ax.set_ylabel("height [m]")
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_field(transects: dict[str, Transect]) -> Figure:
    """The four surveyed depth profiles."""
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, t in transects.items():
        x = t.spacing * np.arange(t.depths.size)
        ax.plot(x, -t.depths, "o-", color=COLORS[name], label=f"{name} (swell {t.swell_direction})")
    ax.axhline(0, color="#0b3d91", lw=1)
    ax.set_xlabel("distance along the transect [m]")
    ax.set_ylabel("sea bed elevation [m]")
    ax.set_title("Field survey: 4 transects, 4 points each")
    ax.legend(frameon=False, ncol=2)
    return fig


def plot_slope_sweep(sweep: pd.DataFrame, band: float) -> Figure:
    """One wave, many plane slopes. ``band`` is the relative numerical uncertainty."""
    fig, ax = plt.subplots(figsize=(7, 4))
    x = sweep["slope_deg"].to_numpy()
    y = sweep["tube_area_norm"].to_numpy(dtype=float)
    ax.fill_between(x, y * (1 - band), y * (1 + band), color="#0072B2", alpha=0.2,
                    label=f"numerical uncertainty ±{band:.0%}")
    ax.plot(x, y, "o-", color="#0072B2", ms=3.5, label="simulated")
    early = sweep["overturned_early"].to_numpy()
    if early.any():
        ax.axvspan(x[early].min() - 0.5, x[early].max() + 0.5, color="#999", alpha=0.15,
                   label="surface folds before the eta/h criterion")
    ax.set_xlabel("slope of the beach [degrees]")
    ax.set_ylabel("tube area / A²  [-]")
    ax.set_title("Same wave (A = 2 m), different slopes")
    ax.legend(frameon=False, fontsize=8)
    return fig


def plot_transects(table: pd.DataFrame) -> Figure:
    """Tube size per transect, and the wave height each transect was given."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    names = table["transect"].tolist()
    colors = [COLORS[n] for n in names]

    ax1.bar(names, table["tube_area_norm"], color=colors)
    for i, row in table.iterrows():
        ax1.text(i, row["tube_area_norm"] + 0.01, f"slope {row['reef_slope_deg']:.2f}°",
                 ha="center", fontsize=8)
    ax1.set_ylabel("tube area / A²  [-]")
    ax1.set_title("Simulated tube size per transect")

    for i, row in table.iterrows():
        ax2.vlines(i, row["window_low"], row["window_high"], color=colors[i], lw=8, alpha=0.35)
        ax2.plot(i, row["amplitude"], "o", color=colors[i], ms=7)
    ax2.set_xticks(range(len(names)), names)
    ax2.set_ylabel("wave height A [m]")
    ax2.set_title("Wave height given to each transect (dot)\nand heights the model accepts (bar)")
    ax2.grid(axis="x", visible=False)
    fig.tight_layout()
    return fig


def plot_width(sens: pd.DataFrame) -> Figure:
    """How the transect comparison changes with the width of the incident wave."""
    fig, ax = plt.subplots(figsize=(7, 4))
    for name, rows in sens.groupby("transect"):
        ax.plot(rows["width"], rows["tube_area_norm"], "o-", color=COLORS[name], label=name)
    ax.set_xlabel("width of the incident wave, sigma [m]")
    ax.set_ylabel("tube area / A²  [-]")
    ax.set_title("The comparison depends on an arbitrary choice")
    ax.legend(frameon=False, ncol=4)
    return fig


def plot_monte_carlo(samples: pd.DataFrame, valid: pd.DataFrame, freq: pd.Series) -> Figure:
    """Left: uncertainty on the slopes. Right: how stable is the ranking?"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    bins = np.linspace(0, 6.5, 40)
    for name, rows in samples.groupby("transect"):
        ax1.hist(rows["reef_slope_deg"], bins=bins, color=COLORS[name], alpha=0.6, label=name)
    ax1.set_xlabel("reef-face slope [degrees]")
    ax1.set_ylabel("number of samples")
    ax1.set_title("Slope uncertainty (rounding errors of 0.1 m)")
    ax1.legend(frameon=False)

    top = freq[freq > 0].head(4)  # only the rankings that actually occurred
    ax2.bar(top.index, top.to_numpy() * 100, color="#555")
    ax2.set_ylabel("share of samples [%]")
    ax2.set_title(
        f"Ranking by tube area / A²\n({valid['sample'].nunique()} usable samples out of "
        f"{samples['sample'].nunique()})"
    )
    return fig


def plot_survey(shares: pd.DataFrame) -> Figure:
    """Preferred swell direction among the 30 surfers."""
    rows = shares[shares["question_id"] == "preferred_swell_direction"]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.barh(rows["option"], rows["count"], color="#0072B2")
    ax.invert_yaxis()
    ax.set_xlabel(f"respondents (out of {int(rows['respondents'].iloc[0])}, several answers allowed)")
    ax.set_title("Preferred swell direction")
    return fig
