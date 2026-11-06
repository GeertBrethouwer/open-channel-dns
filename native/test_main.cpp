// Desktop test for the native solver.
//   c++ -O3 -std=c++17 native/test_main.cpp -o build/ocdns_test            (single precision)
//   c++ -O3 -std=c++17 -DOCDNS_DOUBLE native/test_main.cpp -o build/ocdns_test_d
// Usage: ocdns_test [threads] [T] [Nx Ny Nz]
#include <chrono>
#include <cstdio>
#include <cstdlib>

#include "ocdns.hpp"

using namespace ocdns;

static double fftError(int N) {
  FFT1 f(N);
  std::vector<real> re(N), im(N), tr(N), ti(N), r0(N), i0(N);
  Rng rnd(7);
  for (int i = 0; i < N; ++i) { r0[i] = re[i] = real(rnd() - 0.5); i0[i] = im[i] = real(rnd() - 0.5); }
  // reference DFT in double
  std::vector<double> xr(N), xi(N);
  for (int m = 0; m < N; ++m) {
    double sr = 0, si = 0;
    for (int j = 0; j < N; ++j) {
      const double a = -2 * PI * double(m) * j / N;
      sr += r0[j] * std::cos(a) - i0[j] * std::sin(a);
      si += r0[j] * std::sin(a) + i0[j] * std::cos(a);
    }
    xr[m] = sr; xi[m] = si;
  }
  if (f.pow2) {
    // power-of-two path lives in Solver::fftAxis; check the mixed/DFT line routines only
    return 0;
  }
  f.line(re.data(), im.data(), tr.data(), ti.data(), -1);
  double e = 0;
  for (int m = 0; m < N; ++m) e = std::max(e, std::fabs(re[m] - xr[m]) + std::fabs(im[m] - xi[m]));
  f.line(re.data(), im.data(), tr.data(), ti.data(), +1);
  for (int m = 0; m < N; ++m) e = std::max(e, double(std::fabs(re[m] / N - r0[m]) + std::fabs(im[m] / N - i0[m])));
  return e;
}

int main(int argc, char** argv) {
  const int threads = argc > 1 ? std::atoi(argv[1]) : 0;
  const double T = argc > 2 ? std::atof(argv[2]) : 20;
  Params p;
  if (argc > 5) { p.Nx = std::atoi(argv[3]); p.Ny = std::atoi(argv[4]); p.Nz = std::atoi(argv[5]); }
  p.threads = threads;

  std::printf("precision: %s\n", sizeof(real) == 4 ? "single (float)" : "double");
  for (int N : {6, 10, 12, 15, 20, 24, 30, 36, 40, 45, 48, 60, 75, 80, 96, 120, 7, 14}) {
    FFT1 f(N);
    std::printf("fft N=%-4d %-7s max err %.1e\n", N, f.smooth ? "mixed" : (f.pow2 ? "pow2" : "dft"), fftError(N));
  }

  Solver s(p);
  std::printf("grid %dx%dx%d, threads %d\n", s.Nx, s.Ny, s.Nz, s.threads());
  std::printf("init: div %.2e  Ub %.7f  Re_tau %.1f\n", s.maxDiv(), s.ub(), s.reTau());
  // reference numbers to compare with the JavaScript solver after exactly 20 steps
  for (int i = 0; i < 20; ++i) s.step();
  std::printf("after 20 steps: t %.10f  Re_tau %.10f  tke %.10e\n", s.t, s.reTau(), s.tke());

  auto t0 = std::chrono::steady_clock::now();
  long s0 = s.nstep;
  double next = 0;
  while (s.t < T) {
    s.step();
    if (!s.finite()) { std::printf("NaN at step %ld\n", s.nstep); return 1; }
    if (s.t >= next) {
      std::printf("t=%6.1f dt=%.4f Re_tau=%6.1f tke=%.5f div=%.1e\n", s.t, s.dt, s.reTau(), s.tke(), s.maxDiv());
      next += T / 8;
    }
  }
  const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  std::printf("%.2f ms/step over %ld steps\n", ms / (s.nstep - s0), s.nstep - s0);
  return 0;
}
