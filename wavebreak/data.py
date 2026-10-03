"""Read and check the field data (4 depth transects and a survey of 30 surfers)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

# Depths are rounded to 0.1 m, so "measured - tide = corrected" must hold to this precision.
TIDE_TOLERANCE = 0.051


@dataclass(frozen=True)
class Transect:
    """One surveyed line, listed from offshore (first point) to the reef (last points)."""

    name: str
    swell_direction: str
    depths: np.ndarray  # corrected depth of each point [m]
    spacing: float  # distance between two points [m]

    @property
    def reef_depth(self) -> float:
        """Shallowest depth of the line, taken as the reef crest [m]."""
        return float(self.depths.min())

    @property
    def offshore_depth(self) -> float:
        return float(self.depths[0])


def reef_face_slope_deg(depths: np.ndarray, spacing: float) -> float:
    """Slope of the segment that climbs to the reef crest, in degrees.

    The crest is the shallowest point (never the first one). The slope is measured
    between the crest and the point just offshore of it, as in the original study.
    """
    crest = max(int(np.argmin(depths)), 1)
    rise = depths[crest - 1] - depths[crest]  # metres gained over one spacing
    return float(np.degrees(np.arctan(rise / spacing)))


def check_bathymetry(table: pd.DataFrame) -> None:
    """Stop with a clear message if the survey table is not consistent."""
    problems = []
    if table["point_id"].duplicated().any():
        problems.append("duplicate point_id")
    if (table["corrected_depth_m"] <= 0).any():
        problems.append("a corrected depth is not positive")
    gap = (table["measured_depth_m"] - table["tide_correction_m"] - table["corrected_depth_m"]).abs()
    if (gap > TIDE_TOLERANCE).any():
        problems.append("corrected depth differs from measured depth minus tide")
    if problems:
        raise ValueError("Invalid survey: " + "; ".join(problems))


def load_transects(data_dir: Path = DATA_DIR) -> dict[str, Transect]:
    """Read the survey, check it, and return the transects A-D.

    Depths are recomputed as measured - tide (not copied from the spreadsheet column).
    """
    table = pd.read_csv(data_dir / "bathymetry.csv")
    meta = pd.read_csv(data_dir / "transects.csv").set_index("transect")
    check_bathymetry(table)

    transects = {}
    for name, rows in table.sort_values("point_id").groupby("transect"):
        depths = (rows["measured_depth_m"] - rows["tide_correction_m"]).round(1).to_numpy()
        transects[name] = Transect(
            name=name,
            swell_direction=str(meta.loc[name, "swell_direction"]),
            depths=depths,
            spacing=float(meta.loc[name, "point_spacing_m"]),
        )
    return transects


def load_survey_shares(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Surfer survey: for each answer, how many of the 30 respondents chose it."""
    table = pd.read_csv(data_dir / "survey.csv")
    table["share"] = table["count"] / table["respondents"]
    return table
