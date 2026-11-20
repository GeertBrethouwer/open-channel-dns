# Open Channel DNS

An interactive direct numerical simulation (DNS) of turbulent open channel flow that runs in the browser, on phones as well as desktops. No installation or server-side code is needed.

- **Geometry:** periodic in x (streamwise) and z (spanwise), with a no-slip wall at y = 0 and a free-slip, rigid-lid surface at y = h.
- **Forcing:** constant streamwise mass flux. A uniform streamwise pressure gradient holds U_b = 1; there is no mean spanwise pressure gradient.
- **Settings:** bulk Reynolds number Re_b = U_b h / ν, domain lengths L_x and L_z, and grid sizes N_x, N_y, N_z. Wall stretching, CFL number and the initial disturbances can also be set.
- **Output:** live x–z, x–y and z–y slices of u, u′, v, w and the fluctuating spanwise vorticity ω_z′; mean velocity and Reynolds-stress profiles in wall units; and time histories of Re_τ and the fluctuation energy.

## Numerics

- **Discretisation:** second-order central finite differences on a staggered (MAC) grid, with the convective terms in divergence form. The grid is uniform in x and z and clustered towards the wall in y with a tanh mapping.
- **Time integration:** low-storage third-order Runge–Kutta for convection and the x–z viscous terms, with Crank–Nicolson for the wall-normal viscous term. The default is the four-stage scheme of the SIMSON user guide (Table 5.1) with CFL = 2.5; option `rk: 3` (C++: `Params::rk = 3`) selects the three-stage scheme of Wray. The CFL number is |u|/Δx + |v|/Δy + |w|/Δz.
- **Pressure:** a fractional-step projection at every stage. The Poisson equation is solved exactly with FFTs in x and z (modified wavenumbers) and a tridiagonal solve in y, which keeps the discrete divergence at round-off level. The FFT is radix-2 for powers of two and a self-sorting mixed-radix (Stockham) transform for N = 2^a 3^b 5^c; other sizes fall back to a plain DFT.
- **Initial condition:** a turbulent (1/7 power law) or laminar mean profile, plus large-scale divergence-free modes and weak noise.

The default case is Re_b = 3200 (Re_τ ≈ 200) in a 4h × 2h box on a 64 × 64 × 64 grid (Δx⁺ ≈ 13, Δz⁺ ≈ 6). The resolution is coarse for a second-order scheme, so treat the statistics as indicative. Smaller grids run much faster, especially in the browser.

## Two solvers, one page

| | Website | Android app |
|---|---|---|
| Solver | `src/solver.js`, JavaScript | `native/ocdns.hpp`, C++17 |
| Precision | double | single (`float` fields; sums and statistics in double). Compile with `-DOCDNS_DOUBLE` for double. |
| Threads | 1 | 4 by default, up to the number of cores (setting in the app); work is shared dynamically, so fast and slow phone cores both contribute |
| Runs | in the page's animation loop | on its own thread; the page fetches slices and statistics |

Both solvers implement the same algorithm and start from the same random field. In double precision they agree to at least 10 digits after 20 time steps. The page (`src/app.html`) talks to whichever solver is present: inside the Android app the WebView provides `window.OCDNSNative`, otherwise the JavaScript solver is used. The engine chip in the header shows which one runs.

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
| `native/ocdns.hpp` | C++ port of the solver (header only, multithreaded) |
| `native/test_main.cpp` | Desktop test for the C++ solver (FFT accuracy, divergence, timing) |
| `test/mms.js`, `native/mms_main.cpp` | Manufactured-solution tests of the two solvers (see below) |
| `android/` | Android app: WebView with the page, JNI bridge (`app/src/main/cpp/ocdns_jni.cpp`) to the C++ solver |
| `android/app/src/main/assets/index.html` | The page as used inside the app (written by `build.py`) |

## Building

```bash
python3 build.py
```

To try it locally, serve the folder and open `http://localhost:8000/docs/`:

```bash
python3 -m http.server 8000
```

## C++ solver on the desktop

```bash
mkdir -p build && c++ -O3 -std=c++17 native/test_main.cpp -o build/ocdns_test
```
```bash
./build/ocdns_test 4 100
```

The arguments are the number of threads and the end time; optional Nx Ny Nz follow.

## Verification with a manufactured solution

`test/mms.js` and `native/mms_main.cpp` run both solvers with the same exact solution: three-dimensional, time
dependent and divergence free, with u = v = w = 0 at the wall and v = ∂u/∂y = ∂w/∂y = 0 at the surface (an x–y and
a z–y stream-function mode, an oblique x–z mode, the laminar mean profile and a time-dependent mean spanwise flow).
The body force F = ∂u/∂t + (u·∇)u − ν∇²u is added at every Runge–Kutta stage time through the solver's `forcing`
hook. The exact pressure is zero, so the projection and the flux correction have to remove what F contributes.

```bash
c++ -O3 -std=c++17 -DOCDNS_DOUBLE native/mms_main.cpp -o build/ocdns_mms && ./build/ocdns_mms all 4 1000
```
```bash
node test/mms.js all
```

(Without Node.js: `/System/Library/Frameworks/JavaScriptCore.framework/Versions/Current/Helpers/jsc src/solver.js test/mms.js -- all`.)

Results (Re = 100 for the grid refinement, Re = 1000 for the time refinement; largest error of u, v, w):

| N³ (dt = 10⁻³, t = 0.5) | 8 | 16 | 32 | 64 |
|---|---|---|---|---|
| error | 0.95 | 0.21 | 0.045 | 0.011 |

| dt (16³, t = 1, against dt = 1/2560) | 0.04 | 0.02 | 0.01 | 0.005 | 0.0025 |
|---|---|---|---|---|---|
| 3 stages | 5.6e-4 | 6.6e-5 | 2.2e-5 | 1.1e-5 | 4.9e-6 |
| 4 stages | 7.8e-5 | 3.9e-5 | 1.9e-5 | 9.1e-6 | 4.2e-6 |

- Second order in space; the discrete divergence stays at round-off. The two solvers (C++ in double precision)
  give the same errors to all printed digits, and both RK variants give the same spatial errors.
- In time, both schemes converge at first order at small dt, with the same error. This error is the same at
  Re = 100 and 1000, so it does not come from the viscous splitting. It is most likely the wall slip of the
  non-incremental projection (no pressure gradient in the predictor), which is O(Δt) in both schemes.
- At large dt the third-order Runge–Kutta error dominates. There the four-stage scheme is 7× more accurate.
- `rk: 3` reproduces the earlier solver bitwise (JavaScript and C++, single and double precision).

**Stability** (default case, Re_b = 3200, 64³, C++ in single precision, t = 40–150): both schemes run at CFL 2.0, and the four-stage scheme at the default CFL 2.5.
Above CFL ≈ 3.0 (three stages) or ≈ 3.8 (four stages) Δt stops growing because grid-scale noise raises the maximum
velocity. The three-stage scheme blows up at CFL 5, the four-stage scheme at 8.5. A four-stage step costs about
1.4 times a three-stage step.

## Android app

Requires JDK 17 and the Android SDK with NDK 27.2.12479018 and CMake 3.22.1 (`android/local.properties` points to the SDK). Run `python3 build.py` first so the app gets the current page, then:

```bash
cd android && ./gradlew assembleDebug
```

The APK is written to `android/app/build/outputs/apk/debug/app-debug.apk`.

