// Manufactured-solution (MMS) test of the native solver; the same exact solution as test/mms.js.
//   c++ -O3 -std=c++17 -DOCDNS_DOUBLE native/mms_main.cpp -o build/ocdns_mms
// Usage: ocdns_mms [space|time|all] [threads] [Re]
//
// Exact solution (no-slip wall at y = 0, free-slip surface at y = 1), a = 2 pi/Lx, b = 2 pi/Lz:
//   x-y mode   u =  A1 G'(y) sin(ax),          v = -a A1 G(y) cos(ax)
//   z-y mode   w =  A2 G'(y) sin(bz),          v = -b A2 G(y) cos(bz)
//   x-z mode   u = 2b A3 P(y) cos(ax + 2bz),   w = -a A3 P(y) cos(ax + 2bz)
//   mean       u = 1.5 (2y - y^2),             w = Wm(t) P(y)
// G = y^2 (1 - y)(3 - 2y), P = sin(pi y/2); zero pressure; body force F = du/dt + (u . grad) u - nu lap u
// added at every Runge-Kutta stage time.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include "ocdns.hpp"

using namespace ocdns;

namespace {

struct Term {
  int tf;             // time factor: 0 A1, 1 A2, 2 A3, 3 Wm, 4 one
  double c;           // constant factor
  int yf;             // y factor: 0 G, 1 G', 2 P, 3 U
  bool sn;            // sin (true) or cos
  double kx, kz;
};

struct Exact {
  double a, b;
  std::vector<Term> u, v, w;
  Exact(double Lx, double Lz) : a(2 * PI / Lx), b(2 * PI / Lz) {
    u = {{0, 1, 1, true, a, 0}, {2, 2 * b, 2, false, a, 2 * b}, {4, 1, 3, false, 0, 0}};
    v = {{0, -a, 0, false, a, 0}, {1, -b, 0, false, 0, b}};
    w = {{1, 1, 1, true, 0, b}, {2, -a, 2, false, a, 2 * b}, {3, 1, 2, false, 0, 0}};
  }
  static void T(int f, double t, double& v, double& d) {
    const double e1 = 0.3, e2 = 0.25, e3 = 0.15, e4 = 0.2;
    switch (f) {
      case 0: v = e1 * (1 + 0.5 * std::sin(2 * t)); d = e1 * std::cos(2 * t); break;
      case 1: v = e2 * std::cos(1.5 * t); d = -1.5 * e2 * std::sin(1.5 * t); break;
      case 2: v = e3 * (1 + 0.5 * std::cos(3 * t)); d = -1.5 * e3 * std::sin(3 * t); break;
      case 3: v = e4 * std::sin(t); d = e4 * std::cos(t); break;
      default: v = 1; d = 0;
    }
  }
  static void Y(int f, double y, double r[3]) {
    const double h = PI / 2;
    switch (f) {
      case 0: r[0] = 3 * y * y - 5 * y * y * y + 2 * y * y * y * y; r[1] = 6 * y - 15 * y * y + 8 * y * y * y; r[2] = 6 - 30 * y + 24 * y * y; break;
      case 1: r[0] = 6 * y - 15 * y * y + 8 * y * y * y; r[1] = 6 - 30 * y + 24 * y * y; r[2] = -30 + 48 * y; break;
      case 2: r[0] = std::sin(h * y); r[1] = h * std::cos(h * y); r[2] = -h * h * std::sin(h * y); break;
      default: r[0] = 1.5 * (2 * y - y * y); r[1] = 3 - 3 * y; r[2] = -3;
    }
  }
  // r = value, d/dt, d/dx, d/dy, d/dz, laplacian
  static void eval(const std::vector<Term>& c, double x, double y, double z, double t, double r[6]) {
    for (int q = 0; q < 6; ++q) r[q] = 0;
    for (const Term& m : c) {
      double tv, td, yv[3];
      T(m.tf, t, tv, td); tv *= m.c; td *= m.c;
      Y(m.yf, y, yv);
      const double th = m.kx * x + m.kz * z;
      const double X = m.sn ? std::sin(th) : std::cos(th), dX = m.sn ? std::cos(th) : -std::sin(th);
      const double k2 = m.kx * m.kx + m.kz * m.kz;
      r[0] += tv * yv[0] * X;
      r[1] += td * yv[0] * X;
      r[2] += tv * yv[0] * m.kx * dX;
      r[3] += tv * yv[1] * X;
      r[4] += tv * yv[0] * m.kz * dX;
      r[5] += tv * (yv[2] - k2 * yv[0]) * X;
    }
  }
  double value(const std::vector<Term>& c, double x, double y, double z, double t) const {
    double r[6]; eval(c, x, y, z, t, r); return r[0];
  }
  // component q (0 u, 1 v, 2 w) of the body force
  double force(int q, double nu, double x, double y, double z, double t) const {
    double U[6], V[6], W[6];
    eval(u, x, y, z, t, U); eval(v, x, y, z, t, V); eval(w, x, y, z, t, W);
    const double* c = q == 0 ? U : (q == 1 ? V : W);
    return c[1] + U[0] * c[2] + V[0] * c[3] + W[0] * c[4] - nu * c[5];
  }
};

struct Run {
  std::vector<double> u, v, w;
  double eu, ev, ew, div;
};

Run run(int rk, int N, double Re, double dt, double T, int threads) {
  Params p;
  p.Re = Re; p.Lx = 4; p.Lz = 2; p.Nx = p.Ny = p.Nz = N; p.stretch = 1.5; p.amp = 0; p.rk = rk; p.threads = threads;
  Solver s(p);
  const Exact E(p.Lx, p.Lz);
  const int Nx = s.Nx, Ny = s.Ny, Nz = s.Nz, n = s.n;
  const double dx = s.dx, dz = s.dz, nu = s.nu;
  auto id = [&](int i, int j, int k) { return size_t(i) + size_t(Nx) * k + size_t(n) * j; };
  s.forcing = [&](double t, real* Hu, real* Hv, real* Hw) {
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        for (int i = 0; i < Nx; ++i) {
          Hu[id(i, j, k)] += real(E.force(0, nu, i * dx, s.yc[j], (k + 0.5) * dz, t));
          Hw[id(i, j, k)] += real(E.force(2, nu, (i + 0.5) * dx, s.yc[j], k * dz, t));
          if (j > 0) Hv[id(i, j, k)] += real(E.force(1, nu, (i + 0.5) * dx, s.yf[j], (k + 0.5) * dz, t));
        }
  };
  std::vector<double> u0(size_t(n) * Ny), w0(size_t(n) * Ny), v0(size_t(n) * (Ny + 1), 0.0);
  for (int j = 0; j < Ny; ++j)
    for (int k = 0; k < Nz; ++k)
      for (int i = 0; i < Nx; ++i) {
        u0[id(i, j, k)] = E.value(E.u, i * dx, s.yc[j], (k + 0.5) * dz, 0);
        w0[id(i, j, k)] = E.value(E.w, (i + 0.5) * dx, s.yc[j], k * dz, 0);
        if (j > 0) v0[id(i, j, k)] = E.value(E.v, (i + 0.5) * dx, s.yf[j], (k + 0.5) * dz, 0);
      }
  s.setState(u0, v0, w0, 0.0);
  const long nst = std::lround(T / dt);
  for (long q = 0; q < nst; ++q) s.step(dt);
  Run r;
  s.getState(r.u, r.v, r.w);
  r.eu = r.ev = r.ew = 0;
  for (int j = 0; j < Ny; ++j)
    for (int k = 0; k < Nz; ++k)
      for (int i = 0; i < Nx; ++i) {
        const size_t I = id(i, j, k);
        r.eu = std::max(r.eu, std::fabs(r.u[I] - E.value(E.u, i * dx, s.yc[j], (k + 0.5) * dz, s.t)));
        r.ew = std::max(r.ew, std::fabs(r.w[I] - E.value(E.w, (i + 0.5) * dx, s.yc[j], k * dz, s.t)));
        if (j > 0) r.ev = std::max(r.ev, std::fabs(r.v[I] - E.value(E.v, (i + 0.5) * dx, s.yf[j], (k + 0.5) * dz, s.t)));
      }
  r.div = s.maxDiv();
  return r;
}

double maxDiff(const std::vector<double>& a, const std::vector<double>& b) {
  double m = 0;
  for (size_t I = 0; I < a.size(); ++I) m = std::max(m, std::fabs(a[I] - b[I]));
  return m;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string which = argc > 1 ? argv[1] : "all";
  const int threads = argc > 2 ? std::atoi(argv[2]) : 4;
  const double ReT = argc > 3 ? std::atof(argv[3]) : 1000;
  std::printf("precision: %s\n", sizeof(real) == 4 ? "single (float)" : "double");
  if (which == "space" || which == "all") {
    for (int rk : {3, 4}) {
      std::printf("-- grid refinement, rk = %d, Re = 100, dt = 1e-3, T = 0.5 (max error u, v, w; order)\n", rk);
      double prev = 0;
      for (int N : {8, 16, 32, 64}) {
        const Run r = run(rk, N, 100, 1e-3, 0.5, threads);
        const double m = std::max(r.eu, std::max(r.ev, r.ew));
        std::printf("N = %3d   %.2e  %.2e  %.2e   div %.1e", N, r.eu, r.ev, r.ew, r.div);
        if (prev > 0) std::printf("   order %.2f", std::log(prev / m) / std::log(2.0));
        std::printf("\n");
        prev = m;
      }
    }
  }
  if (which == "time" || which == "all") {
    for (int rk : {3, 4}) {
      const Run ref = run(rk, 16, ReT, 1.0 / 2560, 1, threads);
      std::printf("-- time refinement, rk = %d, Re = %g, 16^3, T = 1: max |u - u_ref| (u, v, w; order)\n", rk, ReT);
      double prev = 0;
      for (double dt : {0.04, 0.02, 0.01, 0.005, 0.0025}) {
        const Run r = run(rk, 16, ReT, dt, 1, threads);
        const double eu = maxDiff(r.u, ref.u), ev = maxDiff(r.v, ref.v), ew = maxDiff(r.w, ref.w);
        const double m = std::max(eu, std::max(ev, ew));
        std::printf("dt = %.4f   %.2e  %.2e  %.2e", dt, eu, ev, ew);
        if (prev > 0) std::printf("   order %.2f", std::log(prev / m) / std::log(2.0));
        std::printf("\n");
        prev = m;
      }
    }
  }
  return 0;
}
