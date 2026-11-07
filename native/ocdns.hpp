// Open channel flow DNS: native C++ port of src/solver.js.
//
// No-slip wall at y = 0, free-slip surface at y = h = 1, periodic in x and z.
// Second-order finite differences on a staggered (MAC) grid, stretched in y.
// Low-storage RK3 (Wray) for convection and x/z diffusion, Crank-Nicolson for
// wall-normal diffusion, fractional-step projection at every stage with an
// FFT (x, z) + tridiagonal (y) Poisson solver. Constant streamwise mass flux
// (U_b = 1); no mean spanwise pressure gradient. Units: h = 1, U_b = 1, nu = 1/Re_b.
//
// Precision: fields and arithmetic in `real` (float by default, double when
// compiled with -DOCDNS_DOUBLE). Sums, statistics and diagnostics use double.
#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace ocdns {

#ifdef OCDNS_DOUBLE
using real = double;
#else
using real = float;
#endif

constexpr double PI = 3.14159265358979323846;

// ---------------------------------------------------------------- thread pool
class Pool {
 public:
  explicit Pool(int n) : n_(std::max(1, n)) {
    for (int t = 1; t < n_; ++t) workers_.emplace_back([this] { loop(); });
  }
  ~Pool() {
    {
      std::lock_guard<std::mutex> lk(m_);
      quit_ = true;
      ++gen_;
    }
    cv_.notify_all();
    for (auto& w : workers_) w.join();
  }
  int size() const { return n_; }
  // Runs f(begin, end) over [0, count) in parallel and returns when all of it is done.
  // Work is shared dynamically: the range is cut into many small chunks that threads
  // take from a shared counter, so fast cores simply process more chunks than slow ones
  // (phones mix fast and slow cores).
  void run(int count, const std::function<void(int, int)>& f) {
    if (n_ == 1 || count < 2) {
      f(0, count);
      return;
    }
    const int chunks = std::min(count, n_ * kChunksPerThread);
    {
      std::lock_guard<std::mutex> lk(m_);
      job_ = &f;
      count_ = count;
      grain_ = (count + chunks - 1) / chunks;
      next_.store(0, std::memory_order_relaxed);
      pending_ = n_ - 1;
      ++gen_;
    }
    cv_.notify_all();
    work(f, count, grain_);
    std::unique_lock<std::mutex> lk(m_);
    done_.wait(lk, [this] { return pending_ == 0; });
    job_ = nullptr;
  }

 private:
  static constexpr int kChunksPerThread = 6;
  void work(const std::function<void(int, int)>& f, int count, int grain) {
    for (;;) {
      const int b = next_.fetch_add(grain, std::memory_order_relaxed);
      if (b >= count) return;
      f(b, std::min(count, b + grain));
    }
  }
  void loop() {
    uint64_t seen = 0;
    for (;;) {
      const std::function<void(int, int)>* job;
      int count, grain;
      {
        std::unique_lock<std::mutex> lk(m_);
        cv_.wait(lk, [&] { return gen_ != seen; });
        seen = gen_;
        if (quit_) return;
        job = job_;
        count = count_;
        grain = grain_;
      }
      work(*job, count, grain);
      {
        std::lock_guard<std::mutex> lk(m_);
        if (--pending_ == 0) done_.notify_one();
      }
    }
  }
  int n_;
  std::vector<std::thread> workers_;
  std::mutex m_;
  std::condition_variable cv_, done_;
  const std::function<void(int, int)>* job_ = nullptr;
  int count_ = 0, pending_ = 0, grain_ = 1;
  std::atomic<int> next_{0};
  uint64_t gen_ = 0;
  bool quit_ = false;
};

// ---------------------------------------------------------------- 1-D FFT
struct FFT1 {
  int N = 1;
  bool pow2 = true;
  std::vector<int> rev;
  std::vector<real> cs, sn;
  explicit FFT1(int n) : N(n), pow2(n > 0 && (n & (n - 1)) == 0) {
    if (pow2) {
      int bits = 0;
      while ((1 << bits) < N) ++bits;
      rev.resize(N);
      for (int i = 0; i < N; ++i) {
        int r = 0;
        for (int b = 0; b < bits; ++b) r |= ((i >> b) & 1) << (bits - 1 - b);
        rev[i] = r;
      }
      cs.resize(std::max(1, N / 2));
      sn.resize(std::max(1, N / 2));
      for (int i = 0; i < N / 2; ++i) {
        cs[i] = real(std::cos(2 * PI * i / N));
        sn[i] = real(std::sin(2 * PI * i / N));
      }
    } else {
      cs.resize(N);
      sn.resize(N);
      for (int i = 0; i < N; ++i) {
        cs[i] = real(std::cos(2 * PI * i / N));
        sn[i] = real(std::sin(2 * PI * i / N));
      }
      // mixed-radix plan when N = 2^a 3^b 5^c
      std::vector<int> radices;
      int r = N;
      for (int p : {5, 3, 2})
        while (r % p == 0) { radices.push_back(p); r /= p; }
      smooth = (r == 1);
      if (smooth) {
        int ncur = N;
        for (int p : radices) {
          Stage st;
          st.p = p;
          st.m = ncur / p;
          st.c.resize(size_t(st.m) * p);
          st.s.resize(size_t(st.m) * p);
          for (int k = 0; k < st.m; ++k)
            for (int t = 0; t < p; ++t) {
              st.c[k * p + t] = real(std::cos(2 * PI * k * t / ncur));
              st.s[k * p + t] = real(std::sin(2 * PI * k * t / ncur));
            }
          for (int j = 0; j < p; ++j) {
            st.rc[j] = real(std::cos(2 * PI * j / p));
            st.rs[j] = real(std::sin(2 * PI * j / p));
          }
          stages.push_back(std::move(st));
          ncur /= p;
        }
      }
    }
  }
  struct Stage {
    int p = 2, m = 1;
    std::vector<real> c, s;
    real rc[5] = {0, 0, 0, 0, 0}, rs[5] = {0, 0, 0, 0, 0};
  };
  bool smooth = false;
  std::vector<Stage> stages;

  // Transform one contiguous line (non-power-of-two N): mixed radix if N = 2^a 3^b 5^c, else DFT.
  // tr/ti are work arrays of length N.
  void line(real* re, real* im, real* tr, real* ti, int sign) const {
    if (smooth) mixed(re, im, tr, ti, sign);
    else dft(re, im, tr, ti, sign);
  }
  // Self-sorting (Stockham) mixed-radix FFT.
  void mixed(real* re, real* im, real* tr, real* ti, int sign) const {
    real *xr = re, *xi = im, *yr = tr, *yi = ti;
    int s = 1;
    real ar[5], ai[5];
    for (const Stage& st : stages) {
      const int p = st.p, m = st.m;
      for (int k = 0; k < m; ++k)
        for (int q = 0; q < s; ++q) {
          for (int r = 0; r < p; ++r) {
            const int o = q + s * (k + r * m);
            ar[r] = xr[o];
            ai[r] = xi[o];
          }
          for (int t = 0; t < p; ++t) {
            real sr = 0, si = 0;
            for (int r = 0; r < p; ++r) {
              const int j = (r * t) % p;
              const real wr = st.rc[j], wi = sign * st.rs[j];
              sr += ar[r] * wr - ai[r] * wi;
              si += ar[r] * wi + ai[r] * wr;
            }
            const real wr = st.c[k * p + t], wi = sign * st.s[k * p + t];
            const int o = q + s * (p * k + t);
            yr[o] = sr * wr - si * wi;
            yi[o] = sr * wi + si * wr;
          }
        }
      std::swap(xr, yr);
      std::swap(xi, yi);
      s *= p;
    }
    if (xr != re) {
      std::copy(xr, xr + N, re);
      std::copy(xi, xi + N, im);
    }
  }
  // Plain DFT of one contiguous line (used when N is not a power of two).
  void dft(real* re, real* im, real* tr, real* ti, int sign) const {
    for (int m = 0; m < N; ++m) {
      real sr = 0, si = 0;
      int idx = 0;
      for (int j = 0; j < N; ++j) {
        const real c = cs[idx], s = sign * sn[idx];
        sr += re[j] * c - im[j] * s;
        si += re[j] * s + im[j] * c;
        idx += m;
        if (idx >= N) idx -= N;
      }
      tr[m] = sr;
      ti[m] = si;
    }
    std::copy(tr, tr + N, re);
    std::copy(ti, ti + N, im);
  }
};

// Mulberry32, identical to the JavaScript version so both start from the same field.
struct Rng {
  uint32_t a;
  explicit Rng(uint32_t s) : a(s) {}
  double operator()() {
    a += 0x6D2B79F5u;
    uint32_t t = a;
    t = (t ^ (t >> 15)) * (t | 1u);
    t ^= t + (t ^ (t >> 7)) * (t | 61u);
    return double(t ^ (t >> 14)) / 4294967296.0;
  }
};

struct Params {
  double Re = 3200, Lx = 4, Lz = 2;
  int Nx = 64, Ny = 64, Nz = 64;
  double stretch = 1.8, cfl = 1.2, amp = 0.3, tavg = 50;
  bool laminar = false;
  uint32_t seed = 1;
  int threads = 0;  // 0 = default (4, or fewer if the device has fewer cores)
};

struct HistoryPoint {
  double t, reTau, tke;
};

struct Profiles {
  bool averaged = false;
  double T = 0, t0 = 0, reTau = 0;
  std::vector<double> yp, Up, Wp, urms, vrms, wrms, uv;
};

// Field codes for slices: 0 = u, 1 = u', 2 = v, 3 = w, 4 = omega_z'
enum Field { F_U = 0, F_UP = 1, F_V = 2, F_W = 3, F_OZ = 4 };

class Solver {
 public:
  explicit Solver(const Params& p) : P(p), pool(pickThreads(p.threads)), fx(p.Nx), fz(p.Nz) {
    Nx = p.Nx; Ny = p.Ny; Nz = p.Nz;
    n = Nx * Nz;
    nu = 1.0 / p.Re;
    Lx = p.Lx; Lz = p.Lz;
    dx = Lx / Nx; dz = Lz / Nz;

    yf.resize(Ny + 1);
    for (int j = 0; j <= Ny; ++j) {
      const double eta = double(j) / Ny;
      yf[j] = p.stretch > 1e-6 ? 1 - std::tanh(p.stretch * (1 - eta)) / std::tanh(p.stretch) : eta;
    }
    yf[0] = 0; yf[Ny] = 1;
    yc.resize(Ny); dyc.resize(Ny); dyf.resize(Ny + 1);
    for (int j = 0; j < Ny; ++j) { yc[j] = 0.5 * (yf[j] + yf[j + 1]); dyc[j] = yf[j + 1] - yf[j]; }
    for (int j = 1; j < Ny; ++j) dyf[j] = yc[j] - yc[j - 1];
    dyf[0] = 2 * yc[0];
    dyf[Ny] = 2 * (1 - yc[Ny - 1]);

    ip.resize(Nx); im.resize(Nx); kp.resize(Nz); km.resize(Nz);
    for (int i = 0; i < Nx; ++i) { ip[i] = (i + 1) % Nx; im[i] = (i - 1 + Nx) % Nx; }
    for (int k = 0; k < Nz; ++k) { kp[k] = (k + 1) % Nz; km[k] = (k - 1 + Nz) % Nz; }

    const size_t N3 = size_t(n) * Ny, N3v = size_t(n) * (Ny + 1);
    u.assign(N3, 0); w.assign(N3, 0); v.assign(N3v, 0);
    Hu.assign(N3, 0); Hw.assign(N3, 0); Hv.assign(N3v, 0);
    Hu1.assign(N3, 0); Hw1.assign(N3, 0); Hv1.assign(N3v, 0);
    du.assign(N3, 0); dw.assign(N3, 0); dv.assign(N3v, 0);
    phi.assign(N3, 0); pr.assign(N3, 0); pim.assign(N3, 0);

    // Poisson: modified wavenumbers, y-operator and its factorisation per (kx, kz)
    std::vector<double> K2(n);
    for (int k = 0; k < Nz; ++k)
      for (int i = 0; i < Nx; ++i) {
        const double kx = (2 / dx) * std::sin(PI * i / Nx), kz = (2 / dz) * std::sin(PI * k / Nz);
        K2[i + Nx * k] = kx * kx + kz * kz;
      }
    pa.resize(Ny); pc.resize(Ny);
    for (int j = 0; j < Ny; ++j) {
      pa[j] = j > 0 ? 1 / (dyc[j] * dyf[j]) : 0;
      pc[j] = j < Ny - 1 ? 1 / (dyc[j] * dyf[j + 1]) : 0;
    }
    cp.assign(N3, 0); pinv.assign(N3, 0);
    {
      std::vector<double> cpd(N3, 0.0);
      for (int j = 0; j < Ny; ++j) {
        const size_t o = size_t(n) * j;
        const double a = pa[j], c = pc[j], b0 = -(pa[j] + pc[j]);
        for (int m = 0; m < n; ++m) {
          if (m == 0 && j == 0) continue;  // phi(k = 0, j = 0) = 0
          const double den = b0 - K2[m] - (j > 0 ? a * cpd[o - n + m] : 0.0);
          pinv[o + m] = real(1 / den);
          cpd[o + m] = c / den;
          cp[o + m] = real(cpd[o + m]);
        }
      }
    }
    mir.resize(n);
    for (int k = 0; k < Nz; ++k)
      for (int i = 0; i < Nx; ++i) mir[i + Nx * k] = ((Nx - i) % Nx) + Nx * ((Nz - k) % Nz);

    dtVisc = 2.0 / (4 * nu * (1 / (dx * dx) + 1 / (dz * dz)));
    planeSum.assign(Ny, 0.0);
    planeMax.assign(Ny + 1, 0.0);
    resetStats();
    avgOn = false;
    initialize();
  }

  // ------------------------------------------------------------ public state
  Params P;
  int Nx, Ny, Nz, n;
  double nu, Lx, Lz, dx, dz;
  std::vector<double> yf, yc, dyc, dyf;
  double t = 0, dt = 0, dpdx = 0, cmax = 0;
  long nstep = 0;
  bool avgOn = false;
  std::vector<HistoryPoint> history;
  int threads() const { return pool.size(); }

  // ------------------------------------------------------------ time step
  void step() {
    const double DT = computeDt(false);
    double dp = 0;
    for (int s = 0; s < 3; ++s) {
      const real g = real(G[s]), z = real(Z[s]), a = real(A[s]), rdt = real(DT);
      computeRHS();
      const size_t N3 = u.size(), nv = size_t(n) * Ny;
      pool.run(int(Ny), [&](int b, int e) {
        for (size_t I = size_t(b) * n; I < size_t(e) * n; ++I) {
          du[I] = rdt * (g * Hu[I] + z * Hu1[I] + a * du[I]);
          dw[I] = rdt * (g * Hw[I] + z * Hw1[I] + a * dw[I]);
        }
        for (int j = std::max(b, 1); j < e; ++j)
          for (size_t I = size_t(j) * n; I < size_t(j + 1) * n; ++I)
            dv[I] = rdt * (g * Hv[I] + z * Hv1[I] + a * dv[I]);
      });
      (void)N3; (void)nv;
      const double beta = 0.5 * A[s] * DT;
      solveImplicitUW(du, beta);
      solveImplicitUW(dw, beta);
      solveImplicitV(dv, beta);
      pool.run(int(Ny), [&](int b, int e) {
        for (size_t I = size_t(b) * n; I < size_t(e) * n; ++I) { u[I] += du[I]; w[I] += dw[I]; }
        for (int j = std::max(b, 1); j < e; ++j)
          for (size_t I = size_t(j) * n; I < size_t(j + 1) * n; ++I) v[I] += dv[I];
      });
      std::swap(Hu, Hu1); std::swap(Hv, Hv1); std::swap(Hw, Hw1);
      project(A[s] * DT);
      dp += A[s] * DT * fixFlux(A[s] * DT);
    }
    dpdx = dp / DT;
    t += DT;
    ++nstep;
    accumulate(DT);
    if (!avgOn && t >= P.tavg) { resetStats(); avgOn = true; }
    if (nstep % 5 == 0) record();
  }

  bool finite() const { return std::isfinite(u[0]) && std::isfinite(v[n]) && std::isfinite(w[0]); }

  // ------------------------------------------------------------ diagnostics
  double planeMean(const std::vector<real>& a, int j) const {
    double s = 0;
    const size_t o = size_t(n) * j;
    for (int m = 0; m < n; ++m) s += a[o + m];
    return s / n;
  }
  double bulk(const std::vector<real>& a) const {
    double ub = 0;
    for (int j = 0; j < Ny; ++j) ub += planeMean(a, j) * dyc[j];
    return ub;
  }
  double ub() const { return bulk(u); }
  double wb() const { return bulk(w); }
  double wallShear() const { return nu * planeMean(u, 0) / yc[0]; }
  double reTau() const { return std::sqrt(std::max(wallShear(), 0.0)) / nu; }
  double maxDiv() const {
    double mx = 0;
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        for (int i = 0; i < Nx; ++i) {
          const size_t I = idx(i, j, k);
          const double d = (u[idx(ip[i], j, k)] - u[I]) / dx + (v[I + n] - v[I]) / dyc[j] +
                           (w[idx(i, j, kp[k])] - w[I]) / dz;
          mx = std::max(mx, std::fabs(d));
        }
    return mx;
  }
  double tke() const {
    double e = 0;
    for (int j = 0; j < Ny; ++j) {
      const size_t o = size_t(n) * j;
      const double mu = planeMean(u, j), mw = planeMean(w, j);
      double s = 0;
      for (int m = 0; m < n; ++m) {
        const double a = u[o + m] - mu, c = w[o + m] - mw, b = 0.5 * (double(v[o + m]) + v[o + n + m]);
        s += a * a + b * b + c * c;
      }
      e += 0.5 * (s / n) * dyc[j];
    }
    return e;
  }

  void resetStats() {
    stT = 0; stT0 = t; stTw = 0;
    stU.assign(Ny, 0); stW.assign(Ny, 0); stuu.assign(Ny, 0); stvv.assign(Ny, 0); stww.assign(Ny, 0); stuv.assign(Ny, 0);
  }
  void startAveragingNow() { resetStats(); avgOn = true; }

  // Averaged statistics once averaging has started, otherwise instantaneous plane averages.
  Profiles profiles() const {
    Profiles out;
    if (avgOn && stT > 0) {
      const double T = stT, ut = std::sqrt(std::max(stTw / T, 1e-30));
      out.averaged = true; out.T = T; out.t0 = stT0; out.reTau = ut / nu;
      for (int j = 0; j < Ny; ++j) {
        out.yp.push_back(yc[j] * ut / nu);
        out.Up.push_back(stU[j] / T / ut);
        out.Wp.push_back(stW[j] / T / ut);
        out.urms.push_back(std::sqrt(std::max(stuu[j] / T, 0.0)) / ut);
        out.vrms.push_back(std::sqrt(std::max(stvv[j] / T, 0.0)) / ut);
        out.wrms.push_back(std::sqrt(std::max(stww[j] / T, 0.0)) / ut);
        out.uv.push_back(-stuv[j] / T / (ut * ut));
      }
      return out;
    }
    const double ut = std::sqrt(std::max(wallShear(), 1e-12));
    out.reTau = ut / nu;
    for (int j = 0; j < Ny; ++j) {
      double su, sv, sw, suu, svv, sww, suv;
      planeMoments(j, su, sv, sw, suu, svv, sww, suv);
      out.yp.push_back(yc[j] * ut / nu);
      out.Up.push_back(su / ut);
      out.Wp.push_back(sw / ut);
      out.urms.push_back(std::sqrt(std::max(suu - su * su, 0.0)) / ut);
      out.vrms.push_back(std::sqrt(std::max(svv - sv * sv, 0.0)) / ut);
      out.wrms.push_back(std::sqrt(std::max(sww - sw * sw, 0.0)) / ut);
      out.uv.push_back(-(suv - su * sv) / (ut * ut));
    }
    return out;
  }

  // ------------------------------------------------------------ slices for display
  // plan: [k*Nx + i] at y index j; side: [j*Nx + i] at z index k; end: [j*Nz + k] at x index i.
  void slices(int field, int si, int sj, int sk, float* plan, float* side, float* end) const {
    std::vector<double> mean(Ny, 0.0);
    if (field == F_UP || field == F_OZ)
      for (int j = 0; j < Ny; ++j) mean[j] = planeMean(u, j);
    for (int k = 0; k < Nz; ++k)
      for (int i = 0; i < Nx; ++i) plan[k * Nx + i] = float(sample(field, i, sj, k, mean));
    for (int j = 0; j < Ny; ++j)
      for (int i = 0; i < Nx; ++i) side[j * Nx + i] = float(sample(field, i, j, sk, mean));
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k) end[j * Nz + k] = float(sample(field, si, j, k, mean));
  }

  // Colour range: [0, max] for u, otherwise +-3 rms over y < 0.6 (same rule as the web app).
  void fieldRange(int field, double& lo, double& hi) const {
    std::vector<double> mean(Ny, 0.0);
    if (field == F_UP || field == F_OZ)
      for (int j = 0; j < Ny; ++j) mean[j] = planeMean(u, j);
    if (field == F_U) {
      double mx = 0;
      for (int j = 0; j < Ny; ++j)
        for (int k = 0; k < Nz; k += 2)
          for (int i = 0; i < Nx; i += 2) mx = std::max(mx, sample(field, i, j, k, mean));
      lo = 0;
      hi = std::max(0.5, std::ceil(mx * 10) / 10);
      return;
    }
    double ss = 0;
    long c = 0;
    for (int j = 0; j < Ny && yc[j] < 0.6; ++j)
      for (int k = 0; k < Nz; k += 2)
        for (int i = 0; i < Nx; i += 2) {
          const double a = sample(field, i, j, k, mean);
          ss += a * a;
          ++c;
        }
    double lim = 3 * std::sqrt(ss / std::max(c, 1L));
    if (!(lim > 0)) lim = 0.1;
    const double mag = std::pow(10.0, std::floor(std::log10(lim)) - 1);  // two significant digits
    lim = std::round(lim / mag) * mag;
    lo = -lim;
    hi = lim;
  }

 private:
  static int pickThreads(int req) {
    if (req > 0) return req;
    const unsigned hc = std::thread::hardware_concurrency();
    return int(std::max(1u, std::min(hc ? hc : 1u, 4u)));
  }
  size_t idx(int i, int j, int k) const { return size_t(i) + size_t(Nx) * k + size_t(n) * j; }

  static constexpr double G[3] = {8.0 / 15, 5.0 / 12, 3.0 / 4};
  static constexpr double Z[3] = {0, -17.0 / 60, -5.0 / 12};
  static constexpr double A[3] = {8.0 / 15, 2.0 / 15, 1.0 / 3};

  Pool pool;
  FFT1 fx, fz;
  std::vector<int> ip, im, kp, km, mir;
  std::vector<real> u, v, w, Hu, Hv, Hw, Hu1, Hv1, Hw1, du, dv, dw, phi, pr, pim, cp, pinv;
  std::vector<double> pa, pc;
  double dtVisc = 1;
  std::vector<double> planeSum, planeMax;
  double stT = 0, stT0 = 0, stTw = 0;
  std::vector<double> stU, stW, stuu, stvv, stww, stuv;

  // ------------------------------------------------------------ initial condition
  void initialize() {
    Rng rnd(P.seed ? P.seed : 1u);
    const double Amp = P.amp;
    auto U = [&](double y) {
      return P.laminar ? 1.5 * (2 * y - y * y) : (8.0 / 7) * std::pow(std::max(y, 0.0), 1.0 / 7);
    };
    const int modes[8][2] = {{0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 1}, {1, 3}, {2, 2}, {1, 0}};
    struct Mode { double a, b, p1, p2, c; };
    auto comp = [&]() {
      std::vector<Mode> ms;
      for (auto& md : modes) {
        Mode m;
        m.a = 2 * PI * md[0] / Lx;
        m.b = 2 * PI * md[1] / Lz;
        m.p1 = 2 * PI * rnd();
        m.p2 = 2 * PI * rnd();
        m.c = (0.5 + rnd()) / std::sqrt(8.0);
        ms.push_back(m);
      }
      return ms;
    };
    const auto mu = comp(), mv = comp(), mw = comp();
    auto field = [](const std::vector<Mode>& ms, double x, double z) {
      double s = 0;
      for (auto& m : ms) s += m.c * std::cos(m.a * x + m.p1) * std::cos(m.b * z + m.p2);
      return s;
    };
    for (int j = 0; j < Ny; ++j) {
      const double env = std::sin(PI * yc[j]);
      for (int k = 0; k < Nz; ++k)
        for (int i = 0; i < Nx; ++i) {
          const size_t I = idx(i, j, k);
          const double xc = (i + 0.5) * dx, zc = (k + 0.5) * dz;
          const double ru = rnd();
          u[I] = real(U(yc[j]) + Amp * env * (field(mu, i * dx, zc) + 0.1 * (ru - 0.5)));
          const double rw = rnd();
          w[I] = real(Amp * env * (field(mw, xc, k * dz) + 0.1 * (rw - 0.5)));
        }
    }
    for (int j = 1; j < Ny; ++j) {
      const double env = std::sin(PI * yf[j]);
      for (int k = 0; k < Nz; ++k)
        for (int i = 0; i < Nx; ++i)
          v[idx(i, j, k)] = real(Amp * env * (field(mv, (i + 0.5) * dx, (k + 0.5) * dz) + 0.1 * (rnd() - 0.5)));
    }
    for (int j = 0; j < Ny; ++j) {
      const double mw0 = planeMean(w, j);
      const size_t o = size_t(n) * j;
      for (int m = 0; m < n; ++m) w[o + m] = real(w[o + m] - mw0);
    }
    project(1);
    fixFlux(1);
    computeDt(true);
    record();
  }

  // ------------------------------------------------------------ right-hand side
  // H = -convection + nu (d2/dx2 + d2/dz2);  L = nu d2/dy2 (stored in du, dv, dw)
  void computeRHS() {
    const real rnu = real(nu), idx_ = real(1 / dx), idz_ = real(1 / dz);
    const real nx2 = real(nu / (dx * dx)), nz2 = real(nu / (dz * dz));
    const real h = real(0.5);
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) {
        const size_t oj = size_t(n) * j, ojp = j < Ny - 1 ? oj + n : oj, ojm = j > 0 ? oj - n : oj;
        const size_t vT = size_t(n) * (j + 1), vB = size_t(n) * j;
        const real ihy = real(1 / dyc[j]);
        const real aLy = real(nu / (dyc[j] * dyf[j])), aUy = real(nu / (dyc[j] * dyf[j + 1]));
        const bool top = j < Ny - 1, bot = j > 0;
        for (int k = 0; k < Nz; ++k) {
          const size_t rk = size_t(Nx) * k, rkP = size_t(Nx) * kp[k], rkM = size_t(Nx) * km[k];
          for (int i = 0; i < Nx; ++i) {
            const size_t iP = ip[i], iM = im[i];
            const size_t I = i + rk + oj;
            // ----- u at (x_i face, yc_j, zc_k)
            const real u0 = u[I], uE = u[iP + rk + oj], uW = u[iM + rk + oj];
            const real uN = u[i + rkP + oj], uS = u[i + rkM + oj];
            const real ucE = h * (u0 + uE), ucW = h * (uW + u0);
            real cu = (ucE * ucE - ucW * ucW) * idx_;
            cu += (h * (u0 + uN) * h * (w[iM + rkP + oj] + w[i + rkP + oj]) -
                   h * (u0 + uS) * h * (w[iM + rk + oj] + w[i + rk + oj])) * idz_;
            const real uT = u[i + rk + ojp], uB = u[i + rk + ojm];
            const real fT = top ? h * (u0 + uT) * h * (v[iM + rk + vT] + v[i + rk + vT]) : real(0);
            const real fB = bot ? h * (u0 + uB) * h * (v[iM + rk + vB] + v[i + rk + vB]) : real(0);
            cu += (fT - fB) * ihy;
            Hu[I] = -cu + nx2 * (uE - 2 * u0 + uW) + nz2 * (uN - 2 * u0 + uS);
            const real ugT = top ? uT : u0, ugB = bot ? uB : -u0;
            du[I] = aUy * (ugT - u0) - aLy * (u0 - ugB);
            // ----- w at (xc_i, yc_j, z_k face)
            const real w0 = w[I], wN = w[i + rkP + oj], wS = w[i + rkM + oj];
            const real wE = w[iP + rk + oj], wW = w[iM + rk + oj];
            const real wcN = h * (w0 + wN), wcS = h * (wS + w0);
            real cw = (wcN * wcN - wcS * wcS) * idz_;
            cw += (h * (u[iP + rkM + oj] + u[iP + rk + oj]) * h * (w0 + wE) -
                   h * (u[i + rkM + oj] + u0) * h * (wW + w0)) * idx_;
            const real wT = w[i + rk + ojp], wB = w[i + rk + ojm];
            const real gT = top ? h * (w0 + wT) * h * (v[i + rkM + vT] + v[i + rk + vT]) : real(0);
            const real gB = bot ? h * (w0 + wB) * h * (v[i + rkM + vB] + v[i + rk + vB]) : real(0);
            cw += (gT - gB) * ihy;
            Hw[I] = -cw + nx2 * (wE - 2 * w0 + wW) + nz2 * (wN - 2 * w0 + wS);
            const real wgT = top ? wT : w0, wgB = bot ? wB : -w0;
            dw[I] = aUy * (wgT - w0) - aLy * (w0 - wgB);
          }
        }
      }
    });
    // ----- v at (xc_i, y_j face, zc_k), j = 1..Ny-1
    pool.run(Ny, [&](int jb, int je) {
      for (int j = std::max(jb, 1); j < je; ++j) {
        const size_t oj = size_t(n) * j, ojm = oj - n;
        const real ihy = real(1 / dyf[j]);
        const real bL = real(nu / (dyf[j] * dyc[j - 1])), bU = real(nu / (dyf[j] * dyc[j]));
        for (int k = 0; k < Nz; ++k) {
          const size_t rk = size_t(Nx) * k, rkP = size_t(Nx) * kp[k], rkM = size_t(Nx) * km[k];
          for (int i = 0; i < Nx; ++i) {
            const size_t iP = ip[i], iM = im[i];
            const size_t I = i + rk + oj;
            const real v0 = v[I], vE = v[iP + rk + oj], vW = v[iM + rk + oj];
            const real vN = v[i + rkP + oj], vS = v[i + rkM + oj];
            const real vT = v[I + n], vB = v[I - n];
            real c = (h * (u[iP + rk + ojm] + u[iP + rk + oj]) * h * (v0 + vE) -
                      h * (u[i + rk + ojm] + u[i + rk + oj]) * h * (vW + v0)) * idx_;
            const real vcT = h * (v0 + vT), vcB = h * (vB + v0);
            c += (vcT * vcT - vcB * vcB) * ihy;
            c += (h * (w[i + rkP + ojm] + w[i + rkP + oj]) * h * (v0 + vN) -
                  h * (w[i + rk + ojm] + w[i + rk + oj]) * h * (vS + v0)) * idz_;
            Hv[I] = -c + nx2 * (vE - 2 * v0 + vW) + nz2 * (vN - 2 * v0 + vS);
            dv[I] = bU * (vT - v0) - bL * (v0 - vB);
          }
        }
      }
    });
    (void)rnu;
  }

  // ------------------------------------------------------------ implicit wall-normal diffusion
  void solveImplicitUW(std::vector<real>& arr, double beta) {
    std::vector<double> lo(Ny, 0), cpv(Ny, 0), den(Ny, 0);
    double prevC = 0;
    for (int j = 0; j < Ny; ++j) {
      const double aL = nu / (dyc[j] * dyf[j]), aU = nu / (dyc[j] * dyf[j + 1]);
      double l = 0, di, up = 0;
      if (Ny == 1) di = 1 + beta * 2 * aL;
      else if (j == 0) { di = 1 + beta * (2 * aL + aU); up = -beta * aU; }
      else if (j == Ny - 1) { di = 1 + beta * aL; l = -beta * aL; }
      else { di = 1 + beta * (aL + aU); l = -beta * aL; up = -beta * aU; }
      den[j] = di - l * prevC;
      cpv[j] = up / den[j];
      prevC = cpv[j];
      lo[j] = l;
    }
    pool.run(n, [&](int mb, int me) {
      for (int j = 0; j < Ny; ++j) {
        const size_t o = size_t(n) * j;
        const real l = real(lo[j]), id = real(1 / den[j]);
        if (j == 0) for (int m = mb; m < me; ++m) arr[m] *= id;
        else for (int m = mb; m < me; ++m) arr[o + m] = (arr[o + m] - l * arr[o - n + m]) * id;
      }
      for (int j = Ny - 2; j >= 0; --j) {
        const size_t o = size_t(n) * j;
        const real c = real(cpv[j]);
        for (int m = mb; m < me; ++m) arr[o + m] -= c * arr[o + n + m];
      }
    });
  }

  void solveImplicitV(std::vector<real>& arr, double beta) {
    if (Ny < 2) return;
    std::vector<double> lo(Ny + 1, 0), cpv(Ny + 1, 0), den(Ny + 1, 1);
    double prevC = 0;
    for (int j = 1; j < Ny; ++j) {
      const double bL = nu / (dyf[j] * dyc[j - 1]), bU = nu / (dyf[j] * dyc[j]);
      const double di = 1 + beta * (bL + bU);
      const double l = j > 1 ? -beta * bL : 0, up = j < Ny - 1 ? -beta * bU : 0;
      den[j] = di - l * prevC;
      cpv[j] = up / den[j];
      prevC = cpv[j];
      lo[j] = l;
    }
    pool.run(n, [&](int mb, int me) {
      for (int j = 1; j < Ny; ++j) {
        const size_t o = size_t(n) * j;
        const real l = real(lo[j]), id = real(1 / den[j]);
        if (j == 1) for (int m = mb; m < me; ++m) arr[o + m] *= id;
        else for (int m = mb; m < me; ++m) arr[o + m] = (arr[o + m] - l * arr[o - n + m]) * id;
      }
      for (int j = Ny - 2; j >= 1; --j) {
        const size_t o = size_t(n) * j;
        const real c = real(cpv[j]);
        for (int m = mb; m < me; ++m) arr[o + m] -= c * arr[o + n + m];
      }
    });
  }

  // ------------------------------------------------------------ FFTs
  // Transform along one axis for `count` lines at once (inner loop over lines).
  void fftAxis(real* re, real* im_, int sign, const FFT1& f, int stride, int count, int cstride,
               std::vector<real>& b1, std::vector<real>& b2, std::vector<real>& b3, std::vector<real>& b4) const {
    const int N = f.N;
    if (N == 1) return;
    if (!f.pow2) {
      b1.resize(N); b2.resize(N); b3.resize(N); b4.resize(N);
      for (int c = 0; c < count; ++c) {
        const size_t o = size_t(c) * cstride;
        for (int a = 0; a < N; ++a) { b1[a] = re[o + size_t(a) * stride]; b2[a] = im_[o + size_t(a) * stride]; }
        f.line(b1.data(), b2.data(), b3.data(), b4.data(), sign);
        for (int a = 0; a < N; ++a) { re[o + size_t(a) * stride] = b1[a]; im_[o + size_t(a) * stride] = b2[a]; }
      }
      return;
    }
    for (int a = 0; a < N; ++a) {
      const int r = f.rev[a];
      if (r > a) {
        const size_t oa = size_t(a) * stride, ob = size_t(r) * stride;
        for (int c = 0; c < count; ++c) {
          const size_t p = oa + size_t(c) * cstride, q = ob + size_t(c) * cstride;
          std::swap(re[p], re[q]);
          std::swap(im_[p], im_[q]);
        }
      }
    }
    for (int size = 2; size <= N; size <<= 1) {
      const int half = size >> 1, stp = N / size;
      for (int k = 0; k < half; ++k) {
        const real wr = f.cs[k * stp], wi = sign * f.sn[k * stp];
        for (int start = 0; start < N; start += size) {
          const size_t oa = size_t(start + k) * stride, ob = oa + size_t(half) * stride;
          for (int c = 0; c < count; ++c) {
            const size_t p = oa + size_t(c) * cstride, q = ob + size_t(c) * cstride;
            const real qr = re[q], qi = im_[q];
            const real xr = qr * wr - qi * wi, xi = qr * wi + qi * wr;
            const real pr_ = re[p], pi_ = im_[p];
            re[q] = pr_ - xr; im_[q] = pi_ - xi;
            re[p] = pr_ + xr; im_[p] = pi_ + xi;
          }
        }
      }
    }
  }
  struct Scratch { std::vector<real> b1, b2, b3, b4, tr, ti; };
  void fft2(real* re, real* im_, int sign, Scratch& s) const {
    fftAxis(re, im_, sign, fx, 1, Nz, Nx, s.b1, s.b2, s.b3, s.b4);   // along x
    fftAxis(re, im_, sign, fz, Nx, Nx, 1, s.b1, s.b2, s.b3, s.b4);   // along z
  }

  // ------------------------------------------------------------ projection
  void project(double adt) {
    const real idx_ = real(1 / dx), idz_ = real(1 / dz), sc = real(1 / adt);
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) {
        const size_t oj = size_t(n) * j;
        const real ihy = real(1 / dyc[j]);
        for (int k = 0; k < Nz; ++k) {
          const size_t rk = size_t(Nx) * k, rkP = size_t(Nx) * kp[k];
          for (int i = 0; i < Nx; ++i) {
            const size_t I = i + rk + oj;
            const real d = (u[ip[i] + rk + oj] - u[I]) * idx_ + (v[I + n] - v[I]) * ihy + (w[i + rkP + oj] - w[I]) * idz_;
            pr[I] = d * sc;
            pim[I] = 0;
          }
        }
      }
    });
    // forward FFTs, two real planes per complex transform
    const int npairs = (Ny + 1) / 2;
    pool.run(npairs, [&](int pb, int pe) {
      Scratch s;
      s.tr.resize(n); s.ti.resize(n);
      for (int q = pb; q < pe; ++q) {
        const int j = 2 * q;
        const size_t o = size_t(n) * j;
        if (j + 1 >= Ny) { fft2(&pr[o], &pim[o], -1, s); continue; }
        const size_t o2 = o + n;
        for (int m = 0; m < n; ++m) pim[o + m] = pr[o2 + m];
        fft2(&pr[o], &pim[o], -1, s);
        for (int m = 0; m < n; ++m) { s.tr[m] = pr[o + m]; s.ti[m] = pim[o + m]; }
        for (int m = 0; m < n; ++m) {
          const real zr = s.tr[m], zi = s.ti[m], mr = s.tr[mir[m]], mi = -s.ti[mir[m]];
          pr[o + m] = real(0.5) * (zr + mr); pim[o + m] = real(0.5) * (zi + mi);
          pr[o2 + m] = real(0.5) * (zi - mi); pim[o2 + m] = real(-0.5) * (zr - mr);
        }
      }
    });
    // tridiagonal solve in y for every (kx, kz)
    pool.run(n, [&](int mb, int me) {
      for (int j = 0; j < Ny; ++j) {
        const size_t o = size_t(n) * j;
        const real a = real(pa[j]);
        if (j == 0) for (int m = mb; m < me; ++m) { pr[m] *= pinv[m]; pim[m] *= pinv[m]; }
        else for (int m = mb; m < me; ++m) {
          const real id = pinv[o + m];
          pr[o + m] = (pr[o + m] - a * pr[o - n + m]) * id;
          pim[o + m] = (pim[o + m] - a * pim[o - n + m]) * id;
        }
      }
      for (int j = Ny - 2; j >= 0; --j) {
        const size_t o = size_t(n) * j;
        for (int m = mb; m < me; ++m) {
          const real c = cp[o + m];
          pr[o + m] -= c * pr[o + n + m];
          pim[o + m] -= c * pim[o + n + m];
        }
      }
    });
    // inverse FFTs, two planes per transform
    const real scl = real(1.0 / n);
    pool.run(npairs, [&](int pb, int pe) {
      Scratch s;
      for (int q = pb; q < pe; ++q) {
        const int j = 2 * q;
        const size_t o = size_t(n) * j;
        if (j + 1 >= Ny) {
          fft2(&pr[o], &pim[o], +1, s);
          for (int m = 0; m < n; ++m) phi[o + m] = pr[o + m] * scl;
          continue;
        }
        const size_t o2 = o + n;
        for (int m = 0; m < n; ++m) {
          const real ar = pr[o + m], ai = pim[o + m], br = pr[o2 + m], bi = pim[o2 + m];
          pr[o + m] = ar - bi;
          pim[o + m] = ai + br;
        }
        fft2(&pr[o], &pim[o], +1, s);
        for (int m = 0; m < n; ++m) { phi[o + m] = pr[o + m] * scl; phi[o2 + m] = pim[o + m] * scl; }
      }
    });
    // correct the velocities
    const real radt = real(adt);
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) {
        const size_t oj = size_t(n) * j;
        for (int k = 0; k < Nz; ++k) {
          const size_t rk = size_t(Nx) * k, rkM = size_t(Nx) * km[k];
          for (int i = 0; i < Nx; ++i) {
            const size_t I = i + rk + oj;
            u[I] -= radt * (phi[I] - phi[im[i] + rk + oj]) * idx_;
            w[I] -= radt * (phi[I] - phi[i + rkM + oj]) * idz_;
          }
        }
        if (j >= 1) {
          const real ihy = real(1 / dyf[j]);
          for (int m = 0; m < n; ++m) v[oj + m] -= radt * (phi[oj + m] - phi[oj - n + m]) * ihy;
        }
      }
    });
  }

  // Constant streamwise mass flux U_b = 1 (uniform mean pressure gradient in x).
  // No mean spanwise pressure gradient. Returns -dP/dx of this substage.
  double fixFlux(double adt) {
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) planeSum[j] = planeMean(u, j) * dyc[j];
    });
    double ub = 0;
    for (int j = 0; j < Ny; ++j) ub += planeSum[j];
    const real cu = real(1 - ub);
    pool.run(Ny, [&](int jb, int je) {
      for (size_t I = size_t(jb) * n; I < size_t(je) * n; ++I) u[I] += cu;
    });
    return (1 - ub) / adt;
  }

  double computeDt(bool init) {
    const double idx_ = 1 / dx, idz_ = 1 / dz;
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) {
        const size_t oj = size_t(n) * j;
        const double ihy = 1 / dyc[j];
        double mx = 1e-12;
        for (int m = 0; m < n; ++m) {
          const size_t I = oj + m;
          const double vv = std::max(std::fabs(double(v[I])), std::fabs(double(v[I + n])));
          const double c = std::fabs(double(u[I])) * idx_ + vv * ihy + std::fabs(double(w[I])) * idz_;
          if (c > mx) mx = c;
        }
        planeMax[j] = mx;
      }
    });
    double mx = 1e-12;
    for (int j = 0; j < Ny; ++j) mx = std::max(mx, planeMax[j]);
    cmax = mx;
    double d = std::min(P.cfl / mx, dtVisc);
    if (!init && dt > 0) d = std::min(d, 1.1 * dt);
    if (init) d *= 0.5;
    dt = d;
    return d;
  }

  void planeMoments(int j, double& su, double& sv, double& sw, double& suu, double& svv, double& sww, double& suv) const {
    const size_t o = size_t(n) * j;
    su = sv = sw = suu = svv = sww = suv = 0;
    for (int m = 0; m < n; ++m) {
      const double a = u[o + m], b = 0.5 * (double(v[o + m]) + v[o + n + m]), c = w[o + m];
      su += a; sv += b; sw += c; suu += a * a; svv += b * b; sww += c * c; suv += a * b;
    }
    su /= n; sv /= n; sw /= n; suu /= n; svv /= n; sww /= n; suv /= n;
  }

  void accumulate(double DT) {
    pool.run(Ny, [&](int jb, int je) {
      for (int j = jb; j < je; ++j) {
        double su, sv, sw, suu, svv, sww, suv;
        planeMoments(j, su, sv, sw, suu, svv, sww, suv);
        stU[j] += DT * su;
        stW[j] += DT * sw;
        stuu[j] += DT * (suu - su * su);
        stvv[j] += DT * (svv - sv * sv);
        stww[j] += DT * (sww - sw * sw);
        stuv[j] += DT * (suv - su * sv);
      }
    });
    stTw += DT * wallShear();
    stT += DT;
  }

  void record() {
    history.push_back({t, reTau(), tke()});
    if (history.size() > 4000) {
      std::vector<HistoryPoint> h;
      for (size_t i = 0; i < history.size(); i += 2) h.push_back(history[i]);
      history.swap(h);
    }
  }

  // Cell-centred value of a display field.
  double sample(int field, int i, int j, int k, const std::vector<double>& mean) const {
    auto uc = [&](int a, int b, int c) { return 0.5 * (double(u[idx(a, b, c)]) + u[idx(ip[a], b, c)]); };
    auto vc = [&](int a, int b, int c) { return 0.5 * (double(v[idx(a, b, c)]) + v[idx(a, b, c) + n]); };
    auto wc = [&](int a, int b, int c) { return 0.5 * (double(w[idx(a, b, c)]) + w[idx(a, b, kp[c])]); };
    switch (field) {
      case F_U: return uc(i, j, k);
      case F_UP: return uc(i, j, k) - mean[j];
      case F_V: return vc(i, j, k);
      case F_W: return wc(i, j, k);
      default: {  // fluctuating spanwise vorticity omega_z' = dv'/dx - du'/dy
        const bool top = j < Ny - 1, bot = j > 0;
        const double uT = top ? uc(i, j + 1, k) - mean[j + 1] : uc(i, j, k) - mean[j];
        const double yT = top ? yc[j + 1] : 2 - yc[j];
        const double uB = bot ? uc(i, j - 1, k) - mean[j - 1] : -(uc(i, j, k) - mean[j]);
        const double yB = bot ? yc[j - 1] : -yc[j];
        return (vc(ip[i], j, k) - vc(im[i], j, k)) / (2 * dx) - (uT - uB) / (yT - yB);
      }
    }
  }
};

}  // namespace ocdns
