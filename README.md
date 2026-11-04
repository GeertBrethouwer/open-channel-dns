# Open Channel DNS

An interactive direct numerical simulation (DNS) of turbulent open channel flow that runs in the browser, on phones as well as desktops. No installation or server-side code is needed.

- **Geometry:** periodic in x (streamwise) and z (spanwise), with a no-slip wall at y = 0 and a free-slip, rigid-lid surface at y = h.
- **Forcing:** constant mass flux. The bulk velocity is held at U_b = 1 by a uniform pressure gradient.
- **Settings:** bulk Reynolds number Re_b = U_b h / ν, domain lengths L_x and L_z, and grid sizes N_x, N_y, N_z. Wall stretching, CFL number and the initial disturbances can also be set.
- **Output:** live x–z, x–y and z–y slices of u, u′, v, w and ω_x; mean velocity and Reynolds-stress profiles in wall units; and time histories of Re_τ and the fluctuation energy.

## Numerics

- **Discretisation:** second-order central finite differences on a staggered (MAC) grid, with the convective terms in divergence form. The grid is uniform in x and z and clustered towards the wall in y with a tanh mapping.
- **Time integration:** three-stage low-storage Runge–Kutta (Wray) for convection and the x–z viscous terms, with Crank–Nicolson for the wall-normal viscous term.
- **Pressure:** a fractional-step projection at every stage. The Poisson equation is solved exactly with FFTs in x and z (modified wavenumbers) and a tridiagonal solve in y, which keeps the discrete divergence at round-off level.
- **Initial condition:** a turbulent (1/7 power law) or laminar mean profile, plus large-scale divergence-free modes and weak noise.

The default case is Re_b = 2000 (Re_τ ≈ 130) in a π h × π/2 h box on a 32 × 48 × 32 grid. It sustains turbulence and runs at roughly 15–20 time steps per second in a desktop browser. The resolution is coarse for a second-order scheme, so treat the statistics as indicative.

## Files

| Path | Purpose |
|---|---|
| `src/solver.js` | The DNS solver (no dependencies) |
| `src/app.html` | Page layout, controls, rendering and plots |
| `build.py` | Inlines the solver into the page and writes the outputs below |
| `docs/index.html` | Standalone page served by GitHub Pages |
| `open-channel-dns.html` | Same page without the document skeleton, for publishing as a Claude artifact |
| `test/test.html` | Browser test harness for the solver (divergence, mass flux, timing) |
| `test/run.js` | The same checks for Node.js |

## Building

```bash
python3 build.py
```

To try it locally, serve the folder and open `http://localhost:8000/docs/`:

```bash
python3 -m http.server 8000
```
