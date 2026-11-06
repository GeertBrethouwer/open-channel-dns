# Open Channel DNS

An interactive direct numerical simulation (DNS) of turbulent open channel flow that runs in the browser, on phones as well as desktops. No installation or server-side code is needed.

- **Geometry:** periodic in x (streamwise) and z (spanwise), with a no-slip wall at y = 0 and a free-slip, rigid-lid surface at y = h.
- **Forcing:** constant streamwise mass flux. A uniform streamwise pressure gradient holds U_b = 1; there is no mean spanwise pressure gradient.
- **Settings:** bulk Reynolds number Re_b = U_b h / ν, domain lengths L_x and L_z, and grid sizes N_x, N_y, N_z. Wall stretching, CFL number and the initial disturbances can also be set.
- **Output:** live x–z, x–y and z–y slices of u, u′, v, w and the fluctuating spanwise vorticity ω_z′; mean velocity and Reynolds-stress profiles in wall units; and time histories of Re_τ and the fluctuation energy.

## Numerics

- **Discretisation:** second-order central finite differences on a staggered (MAC) grid, with the convective terms in divergence form. The grid is uniform in x and z and clustered towards the wall in y with a tanh mapping.
- **Time integration:** three-stage low-storage Runge–Kutta (Wray) for convection and the x–z viscous terms, with Crank–Nicolson for the wall-normal viscous term.
- **Pressure:** a fractional-step projection at every stage. The Poisson equation is solved exactly with FFTs in x and z (modified wavenumbers) and a tridiagonal solve in y, which keeps the discrete divergence at round-off level. The FFT is radix-2 for powers of two and a self-sorting mixed-radix (Stockham) transform for N = 2^a 3^b 5^c; other sizes fall back to a plain DFT.
- **Initial condition:** a turbulent (1/7 power law) or laminar mean profile, plus large-scale divergence-free modes and weak noise.

The default case is Re_b = 2000 (Re_τ ≈ 130) in a π h × π/2 h box on a 32 × 48 × 32 grid. It sustains turbulence and runs at roughly 15–20 time steps per second in a desktop browser. The resolution is coarse for a second-order scheme, so treat the statistics as indicative.

## Two solvers, one page

| | Website | Android app |
|---|---|---|
| Solver | `src/solver.js`, JavaScript | `native/ocdns.hpp`, C++17 |
| Precision | double | single (`float` fields; sums and statistics in double). Compile with `-DOCDNS_DOUBLE` for double. |
| Threads | 1 | up to 4 |
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

## Android app

Requires JDK 17 and the Android SDK with NDK 27.2.12479018 and CMake 3.22.1 (`android/local.properties` points to the SDK). Run `python3 build.py` first so the app gets the current page, then:

```bash
cd android && ./gradlew assembleDebug
```

The APK is written to `android/app/build/outputs/apk/debug/app-debug.apk`.

