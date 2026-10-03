# How the model works

This page explains the model step by step. Everything described here is in
[`wavebreak/model.py`](../wavebreak/model.py) (about 270 lines, comments included).

## 1. The question

A surfer's "tube" (barrel) appears when the top of a wave, the *lip*, is thrown
forward and lands in front of the wave. The model asks: **given the shape of the
sea bed, does the wave break this way, and how big is the tube?**

## 2. The ingredients

| Ingredient | In the code | Meaning |
|---|---|---|
| Sea bed | `depth(x)` | water depth `h` at position `x` (m). Either a plane beach or the 4 surveyed profiles (linear interpolation between the points). |
| Wave | `x`, `eta` arrays | the water surface is described by 900 points ("markers"). Marker `i` is at horizontal position `x[i]` and height `eta[i]` above still water. |
| Start | a Gaussian hump | `eta = A · exp(−(x − x_crest)² / (2σ²))`: height `A`, width `σ`, placed offshore of the sea bed. |

Positions are in metres, times in seconds, `g = 9.81 m/s²`.

## 3. What happens at each time step (dt = 0.01 s)

**Step 1: propagation.** Every marker moves toward the shore at the shallow-water
wave speed

```
c = sqrt( g · (h + eta) )        x ← x + c · dt
```

* where it is deep, the wave goes faster; where it is shallow, slower;
* a marker with a large `eta` (the crest) goes faster than a low one: **the crest
  catches up with the front of the wave, which therefore gets steeper**;
* the `eta` term is only used from `x = 0` onward (offshore the wave is treated as
  a linear wave, as in the original script).

**Step 2: has the wave reached the breaking limit?** The wave starts to break when
its height becomes a large fraction of the local depth:

```
eta_max / h ≥ γ = 0.78
```

This is a *breaker index* criterion. The value 0.78 comes from solitary-wave theory;
measured values range from about 0.6 to 0.8 for gentle (spilling) breakers and 0.8 to
1.2 for plunging ones, and increase with the slope of the bed [1]. The original script
cites Dean & Dalrymple (1991) for 0.78 (not checked again here). The first time the
criterion is true is the **onset** of breaking.

**Step 3: release of the lip.** After onset, a marker that has caught up with the
marker in front of it, and is high enough (more than 4 % of `A`), would make the
surface fold over itself. It is "released": from now on it flies freely, with

```
horizontal speed = (1 + b) · c       with b = 0.4
vertical speed   = 0 at release
```

The factor `b = 0.4` makes the jet faster than the wave. In the original script it stands
for the orbital velocity of the water at the crest (typical values tried: 0 to 0.6). It is
**empirical**: it is not derived from a law (see limitations).

**Step 4: free fall.** A released marker is a projectile:

```
vz ← vz − g · dt      eta ← eta + vz · dt      x ← x + vx · dt
```

The first marker that reaches `eta ≤ 0` (the still-water level) closes the tube.
The run stops a few steps later.

## 4. What is measured

```
        lip (released at x_r, height H)
          ╲
   onset   ╲___      L = distance flown by the lip
     x_o       ╲___
                   ╲ x_c  (closure: lip touches the water)
```

| Quantity | Definition | Comment |
|---|---|---|
| `jet_reach` | `x_c − x_o` | the metric of the original TIPE |
| `lip_throw` | `x_c − x_r` (= `L`) | distance flown by the closing lip |
| `tube_area` | `2/3 · H · L` | area of air under the parabola of the lip (m²) |
| `tube_area_norm` | `tube_area / A²` | no unit; compares waves of different sizes |

Why `2/3 · H · L`: a point launched horizontally from height `H` follows a parabola;
the area between the parabola and the water is exactly two thirds of the rectangle
`H × L`.

## 5. Parameters and where they come from

| Parameter | Value | Status |
|---|---|---|
| `γ` (breaking limit) | 0.78 | literature (solitary-wave value), **fixed**; real values vary from about 0.6 to 1.2 with the slope and the type of breaker [1] |
| `b` (jet boost) | 0.4 | **empirical**, not derived |
| `σ` (width of the wave) | 3 m (field), 4 m (beach) | **a choice**; the results depend on it a lot (see README) |
| `A` (height of the wave) | `1.02 · γ · h_reef` (field), 2 m (beach) | field: rule of the original study; beach: arbitrary |
| `dt`, number of markers | 0.01 s, 900 (field) / 800 (beach) | numerical; changing them moves the tube area by a few % |

## 6. When is the model usable? (validity window)

If the initial wave is already higher than `γ · h` where it is placed, it "breaks"
before the experiment starts, and the result is meaningless. For a transect the model
is usable only for

```
γ · h_reef ≤ A < γ · h_offshore
```

(lower bound: the crest must be able to reach the criterion above the reef).
`starts_breaking(...)` in the code tests the upper bound; the Monte Carlo discards the
samples where it fails.

## 7. How it was checked

| Test | What it shows |
|---|---|
| Speed of a tiny wave on a flat bed = `sqrt(g h)` | propagation step is correct |
| Time at which the front overturns on a flat bed vs. a formula (method of characteristics) | the steepening is correct: 0.04 % difference with the default grid |
| Distance flown by the lip vs. `u · sqrt(2H/g)` | the free fall is correct: 0.4 % difference |
| `2/3 · H · L` on a hand-made example | the area formula |
| Same numbers as the original script logic and the first version of this project | the simplification did not change the results |

All these tests are in [`tests/`](../tests/). They check that the code solves **its own
equations** correctly; they do **not** prove that the equations describe real waves.
There is no comparison with a measured breaking wave.

## Reference

[1] J. Bosboom and M. J. F. Stive, *Coastal Dynamics*, Delft University of Technology
(open textbook), section 5.2.5 "Wave breaking":
<https://geo.libretexts.org/Bookshelves/Oceanography/Coastal_Dynamics_(Bosboom_and_Stive)/05%3A_Coastal_hydrodynamics/5.02%3A_Wave_transformation/5.2.5%3A_Wave_breaking>
(consulted on 2026-10-03).
