# Limitations

Most important first. A result of this project should always be read with these in mind.

## About the model

1. **No validation against a real wave.** The tests show that the code does what the
   equations say. Nothing here shows that the equations reproduce a real breaking wave.
2. **The wave width `σ` is a choice and it matters.** With `σ = 3 m` the ranking of the
   transects is C > D > B > A; with `σ = 4 m` it is D > C > B > A; with `σ = 5 m` it is back to
   C > D > B > A; with `σ = 6 m`, A > C > D > B (README, section 3). A real swell is a train of long waves, not one 3 m bump.
3. **Two criteria for breaking disagree.** The model uses (a) the height / depth ratio
   (`eta/h ≥ 0.78`) and (b) the surface folding over. On gentle slopes (≤ 6° on the plane
   beach, and on all four transects) the surface folds *before* (a) is satisfied. In that
   regime the metrics are not reliable; `jet_reach` can even be negative.
4. **Empirical jet speed.** `b = 0.4` is not derived from physics, and the tube area
   depends on it.
5. **The wave does not grow when the water gets shallower** (no shoaling): `eta` stays
   constant on each marker. The wave height at breaking is imposed, not predicted.
6. **The lip lands on still water (`eta = 0`)**, not on the moving surface in front of it.
7. **2D only, one direction, no wind, current or friction.** The four transects are
   treated as independent problems. Real waves refract around a reef.
8. **Breaking limit fixed at 0.78.** It was not varied. Textbooks give 0.8–1.2 for plunging breakers and a value that increases with the slope (see the reference in [`model.md`](model.md)), so 0.78 probably makes the model break too early on steep slopes.

## About the field data

1. **Only 4 points per transect**, depths rounded to 0.1 m. One reading step over 12 m is
   0.48°: transects B and C have exactly the same slope, and D's slope is one step.
2. **Spacing is "about 12 m"**, the real positions were not recorded. It is not part of the
   Monte Carlo (it was in a previous, more complex version of this project).
3. **The tide correction column** is in metres even though it was called a coefficient in
   the spreadsheet; it is interpreted as a water-level correction because
   `corrected = measured − correction` holds for all 16 points.
4. **A single survey**: no repeat measurement, no date or instrument recorded.

## About the comparison between transects

* To compare slopes fairly the wave must be the same on every transect. It cannot be:
  no wave height is accepted by the model on all four transects (README, figure
  "transects", right). The original study gave each transect its own wave height
  (`1.02 · γ · h_reef`), so **deeper reefs get bigger waves**: the comparison mixes the
  effect of slope with the effect of wave size.
* The only fair test of the slope is done on an idealised plane beach (README, "One wave,
  many slopes").

## About the numerics

* Changing the number of markers and the time step moves the tube area by about ±4 %
  (3 resolutions tested, one case: 10° beach, README). Smaller differences are not
  meaningful.
* The metrics depend on which discrete marker closes the tube: the sweeps have small jumps.

## About the surfer survey

* 30 people, convenience sample, only aggregated answers (no individual answers).
* People were asked their **preferred swell direction**, which is not a measurement of the
  *shape of a wave*; it also reflects how often a swell occurs, its size, the surfer's level.
* Four directions cannot give a statistically significant ranking, so the survey is shown
  but **not compared quantitatively** with the model.

## What this project does and does not show

**It shows**: how a simple kinematic model of a breaking wave can be built, tested against
formulas, and used for numerical experiments; that, in this model, slope has a
**non-monotonic** effect on the tube; and that the first ranking of the transects is not
robust.

**It does not show**: which of the four surveyed reefs makes the best waves, nor any
quantitative prediction for a real surf spot.
