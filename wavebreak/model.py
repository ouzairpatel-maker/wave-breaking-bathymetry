"""Surface-marker model of a wave that steepens and breaks (plunging breaker).

The idea in four steps (the physics is explained in docs/model.md):

1. The water surface is a line of N points ("markers"), each with a position x
   and an elevation eta above the still-water level.
2. Every marker moves toward the shore at the shallow-water wave speed
   c = sqrt(g * (h + eta)), where h is the local depth. A high marker (crest)
   goes faster than a low one, so the front of the wave gets steeper.
3. Once the wave is high enough compared with the depth (eta / h >= gamma),
   a marker that catches up with the marker in front of it is released: the
   surface would fold over, which is the lip of the breaking wave.
4. A released marker flies like a projectile (constant horizontal speed, free
   fall). The first one to touch the water closes the "tube" (the barrel).

Elevation eta is measured from the still-water level (eta = 0), x in metres,
time in seconds.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

import numpy as np

G = 9.81  # gravity, m/s^2

# State of a marker
PROPAGATING, FALLING, LANDED = 0, 1, 2

# Lip contacts during the first steps after onset do not close the tube, and the
# run goes on a few steps after closure so that the last frames can be drawn.
CLOSURE_GUARD_STEPS = 5
POST_CLOSURE_STEPS = 5

DepthFunction = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Params:
    """Numbers that control a simulation. Defaults = the field study of the TIPE."""

    width: float = 3.0  # standard deviation sigma of the initial Gaussian wave [m]
    crest_offset: float = 15.0  # initial crest position, this far offshore of x = 0 [m]
    gamma: float = 0.78  # breaking threshold eta / h (McCowan) [-]
    jet_boost: float = 0.4  # extra jet speed at release, fraction of c (empirical) [-]
    lip_threshold: float | None = None  # a marker must be this high to be released [m]
    nonlinear_from: float | None = 0.0  # x where the crest starts to go faster (None = everywhere)
    dt: float = 0.01  # time step [s]
    n_markers: int = 900
    x_min: float = -45.0
    x_max: float = 110.0
    max_steps: int = 12000


FIELD = Params()
BEACH = Params(
    width=4.0,
    crest_offset=20.0,
    lip_threshold=0.2,
    n_markers=800,
    x_min=-55.0,
    x_max=130.0,
)


# --------------------------------------------------------------------- sea bed


def plane_beach(offshore_depth: float, slope_deg: float, toe_x: float = 0.0) -> DepthFunction:
    """Flat bed of depth ``offshore_depth``, then a plane slope starting at ``toe_x``."""
    tan_beta = np.tan(np.radians(slope_deg))

    def depth(x: np.ndarray) -> np.ndarray:
        return offshore_depth - tan_beta * np.maximum(np.asarray(x) - toe_x, 0.0)

    return depth


def survey_profile(depths: np.ndarray, spacing: float) -> DepthFunction:
    """Depth linearly interpolated between survey points, flat outside the surveyed range."""
    x_points = spacing * np.arange(len(depths))

    def depth(x: np.ndarray) -> np.ndarray:
        return np.interp(np.asarray(x), x_points, depths)

    return depth


# ---------------------------------------------------------------------- result


@dataclass
class Result:
    """What happened during one run. Steps are counted from 0; time = (step + 1) * dt."""

    amplitude: float
    dt: float
    onset_step: int | None = None  # eta / h >= gamma reached
    onset_x: float | None = None
    closure_step: int | None = None  # first lip contact with the water
    closure_x: float | None = None
    release_x: float | None = None  # where the closing lip was released
    release_height: float | None = None  # its height above still water
    first_crossing_step: int | None = None  # first time the surface folds over
    x_history: list[np.ndarray] = field(default_factory=list)  # filled if record=True
    eta_history: list[np.ndarray] = field(default_factory=list)

    @property
    def broke(self) -> bool:
        """True if a tube closed (the lip reached the water)."""
        return self.closure_step is not None

    @property
    def jet_reach(self) -> float | None:
        """Original TIPE metric: closure position minus onset position [m]."""
        if self.closure_x is None or self.onset_x is None:
            return None
        return self.closure_x - self.onset_x

    @property
    def lip_throw(self) -> float | None:
        """Horizontal distance flown by the closing lip, from release to landing [m]."""
        if self.closure_x is None or self.release_x is None:
            return None
        return self.closure_x - self.release_x

    @property
    def tube_area(self) -> float | None:
        """Air area under the lip's flight, 2/3 * height * throw [m^2].

        A projectile launched horizontally from height H lands a distance L away.
        The area between its parabola and the water is exactly 2/3 * H * L.
        """
        if self.lip_throw is None or self.release_height is None:
            return None
        return 2.0 / 3.0 * self.release_height * self.lip_throw

    @property
    def tube_area_normalized(self) -> float | None:
        """Tube area divided by amplitude squared, to compare waves of different size [-]."""
        area = self.tube_area
        return None if area is None else area / self.amplitude**2

    @property
    def overturned_early(self) -> bool:
        """True if the surface folded over BEFORE the eta / h >= gamma criterion was met.

        The model has two breaking criteria that do not always agree. This flag shows
        when they disagree (typical on gentle slopes).
        """
        if self.first_crossing_step is None:
            return False
        return self.onset_step is None or self.first_crossing_step < self.onset_step

    def time(self, step: int) -> float:
        """Physical time of the state computed at ``step`` [s]."""
        return (step + 1) * self.dt


def starts_breaking(depth: DepthFunction, amplitude: float, params: Params) -> bool:
    """True if the initial wave already satisfies eta / h >= gamma where it is placed.

    The model is not valid then (the wave breaks before the experiment starts).
    """
    h = float(depth(np.array(-params.crest_offset)))
    return amplitude >= params.gamma * h


# ------------------------------------------------------------------ simulation


def simulate(
    depth: DepthFunction, amplitude: float, params: Params = FIELD, record: bool = False
) -> Result:
    """Run the model with a Gaussian incident wave of crest height ``amplitude`` [m].

    Set ``record=True`` to keep every frame (needed to draw the wave, slower).
    """
    p = params
    dt = p.dt
    lip_min = p.lip_threshold if p.lip_threshold is not None else 0.04 * amplitude

    # Initial state: Gaussian hump on still water.
    x = np.linspace(p.x_min, p.x_max, p.n_markers)
    eta = amplitude * np.exp(-((x + p.crest_offset) ** 2) / (2.0 * p.width**2))
    vx = np.zeros_like(x)
    vz = np.zeros_like(x)
    state = np.full(x.size, PROPAGATING)
    release_x = np.full(x.size, np.nan)
    release_eta = np.full(x.size, np.nan)

    res = Result(amplitude=amplitude, dt=dt)
    if record:
        res.x_history.append(x.copy())
        res.eta_history.append(eta.copy())

    for step in range(p.max_steps):
        # 1) Propagation: wave speed from the local depth (and the height, once nonlinear).
        moving = state == PROPAGATING
        if p.nonlinear_from is None:
            uses_eta = np.ones_like(x)
        else:
            uses_eta = (x >= p.nonlinear_from).astype(float)
        c = np.sqrt(G * np.maximum(depth(x) + uses_eta * eta, 0.0))
        x[moving] += dt * c[moving]

        # Has the surface folded over (a marker passed the one in front of it)?
        if res.first_crossing_step is None:
            folded = (x[:-1] >= x[1:]) & (eta[:-1] > lip_min)
            if np.any(folded):
                res.first_crossing_step = step

        # 2) Breaking onset: highest marker compared with the local depth.
        if res.onset_step is None:
            top = int(np.argmax(eta))
            h_top = float(depth(x[top]))
            if h_top > 0.0 and eta[top] >= p.gamma * h_top:
                res.onset_step, res.onset_x = step, float(x[top])

        # 3) Lip release: after onset, a marker that caught up with its neighbour is thrown.
        if res.onset_step is not None:
            active = np.flatnonzero(state == PROPAGATING)
            if active.size >= 2:
                caught_up = np.zeros(active.size, dtype=bool)
                caught_up[:-1] = x[active][:-1] >= x[active][1:]
                caught_up &= eta[active] > lip_min
                released = active[caught_up]
                if released.size:
                    state[released] = FALLING
                    vx[released] = (1.0 + p.jet_boost) * c[released]
                    vz[released] = 0.0
                    release_x[released] = x[released]
                    release_eta[released] = eta[released]

        # 4) Free fall of the released markers (Euler step), landing on still water.
        falling = state == FALLING
        if np.any(falling):
            vz[falling] -= G * dt
            eta[falling] += vz[falling] * dt
            x[falling] += vx[falling] * dt
            touching = falling & (eta <= 0.0)
            first_contact = (
                res.onset_step is not None
                and res.closure_step is None
                and step > res.onset_step + CLOSURE_GUARD_STEPS
                and np.any(touching)
            )
            if first_contact:
                candidates = np.flatnonzero(touching)
                closing = int(candidates[np.argmax(x[candidates])])  # most shoreward one
                res.closure_step, res.closure_x = step, float(x[closing])
                res.release_x = float(release_x[closing])
                res.release_height = float(release_eta[closing])
            state[touching] = LANDED

        if record:
            res.x_history.append(x.copy())
            res.eta_history.append(eta.copy())

        if res.closure_step is not None and step > res.closure_step + POST_CLOSURE_STEPS:
            break

    return res
