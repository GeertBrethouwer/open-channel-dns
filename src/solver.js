"use strict";
/*
 * Open channel flow DNS: no-slip wall at y=0, free-slip surface at y=h=1.
 * Periodic in x (streamwise) and z (spanwise).
 * Second-order finite differences on a staggered (MAC) grid, stretched in y.
 * Time integration: low-storage RK3 (Wray) for convection and x/z diffusion,
 * Crank-Nicolson for wall-normal diffusion, fractional-step projection per substage.
 * Pressure Poisson: FFT in x and z (modified wavenumbers) + tridiagonal solve in y.
 * Constant mass flux: bulk velocity U_b = 1 enforced after every substage.
 * Units: h = 1, U_b = 1, nu = 1/Re_b.
 */

class FFT1 {
  constructor(N) {
    this.N = N;
    this.pow2 = N > 0 && (N & (N - 1)) === 0;
    this.re = new Float64Array(N);
    this.im = new Float64Array(N);
    this.tr = new Float64Array(N);
    this.ti = new Float64Array(N);
    if (this.pow2) {
      const bits = Math.round(Math.log2(N));
      this.rev = new Uint32Array(N);
      for (let i = 0; i < N; i++) {
        let r = 0;
        for (let b = 0; b < bits; b++) r |= ((i >> b) & 1) << (bits - 1 - b);
        this.rev[i] = r;
      }
      this.cs = new Float64Array(Math.max(1, N >> 1));
      this.sn = new Float64Array(Math.max(1, N >> 1));
      for (let i = 0; i < N >> 1; i++) {
        this.cs[i] = Math.cos((2 * Math.PI * i) / N);
        this.sn[i] = Math.sin((2 * Math.PI * i) / N);
      }
    } else {
      this.cs = new Float64Array(N);
      this.sn = new Float64Array(N);
      for (let i = 0; i < N; i++) {
        this.cs[i] = Math.cos((2 * Math.PI * i) / N);
        this.sn[i] = Math.sin((2 * Math.PI * i) / N);
      }
      // mixed-radix plan when N = 2^a 3^b 5^c
      const radices = [];
      let r = N;
      for (const p of [5, 3, 2]) while (r % p === 0) { radices.push(p); r /= p; }
      this.smooth = r === 1;
      if (this.smooth) {
        this.stages = [];
        let ncur = N;
        for (const p of radices) {
          const m = ncur / p, c = new Float64Array(m * p), s = new Float64Array(m * p);
          for (let k = 0; k < m; k++) for (let t = 0; t < p; t++) {
            c[k * p + t] = Math.cos((2 * Math.PI * k * t) / ncur);
            s[k * p + t] = Math.sin((2 * Math.PI * k * t) / ncur);
          }
          const rc = new Float64Array(p), rs = new Float64Array(p);
          for (let j = 0; j < p; j++) { rc[j] = Math.cos((2 * Math.PI * j) / p); rs[j] = Math.sin((2 * Math.PI * j) / p); }
          this.stages.push({ p, m, c, s, rc, rs });
          ncur = m;
        }
        this.ar = new Float64Array(5); this.ai = new Float64Array(5);
      }
    }
  }
  // Self-sorting (Stockham) mixed-radix FFT for N = 2^a 3^b 5^c, using tr/ti as work arrays.
  mixed(sign) {
    let xr = this.re, xi = this.im, yr = this.tr, yi = this.ti, s = 1;
    const ar = this.ar, ai = this.ai;
    for (const st of this.stages) {
      const { p, m, c, s: sn, rc, rs } = st;
      for (let k = 0; k < m; k++) {
        for (let q = 0; q < s; q++) {
          for (let r = 0; r < p; r++) { const o = q + s * (k + r * m); ar[r] = xr[o]; ai[r] = xi[o]; }
          for (let t = 0; t < p; t++) {
            let sr = 0, si = 0;
            for (let r = 0; r < p; r++) {
              const j = (r * t) % p, wr = rc[j], wi = sign * rs[j];
              sr += ar[r] * wr - ai[r] * wi; si += ar[r] * wi + ai[r] * wr;
            }
            const wr = c[k * p + t], wi = sign * sn[k * p + t], o = q + s * (p * k + t);
            yr[o] = sr * wr - si * wi; yi[o] = sr * wi + si * wr;
          }
        }
      }
      let tmp = xr; xr = yr; yr = tmp; tmp = xi; xi = yi; yi = tmp;
      s *= p;
    }
    if (xr !== this.re) { this.re.set(xr); this.im.set(xi); }
  }
  // In-place transform of this.re/this.im. sign = -1 forward, +1 inverse (unscaled).
  run(sign) {
    const N = this.N, re = this.re, im = this.im;
    if (this.pow2) {
      const rev = this.rev;
      for (let i = 0; i < N; i++) {
        const r = rev[i];
        if (r > i) {
          let t = re[i]; re[i] = re[r]; re[r] = t;
          t = im[i]; im[i] = im[r]; im[r] = t;
        }
      }
      for (let size = 2; size <= N; size <<= 1) {
        const half = size >> 1, step = N / size;
        for (let start = 0; start < N; start += size) {
          for (let k = 0; k < half; k++) {
            const wr = this.cs[k * step], wi = sign * this.sn[k * step];
            const a = start + k, b = a + half;
            const xr = re[b] * wr - im[b] * wi;
            const xi = re[b] * wi + im[b] * wr;
            re[b] = re[a] - xr; im[b] = im[a] - xi;
            re[a] += xr; im[a] += xi;
          }
        }
      }
    } else if (this.smooth) {
      this.mixed(sign);
    } else {
      const tr = this.tr, ti = this.ti, cs = this.cs, sn = this.sn;
      for (let m = 0; m < N; m++) {
        let sr = 0, si = 0, idx = 0;
        for (let j = 0; j < N; j++) {
          const c = cs[idx], s = sign * sn[idx];
          sr += re[j] * c - im[j] * s;
          si += re[j] * s + im[j] * c;
          idx += m; if (idx >= N) idx -= N;
        }
        tr[m] = sr; ti[m] = si;
      }
      re.set(tr); im.set(ti);
    }
  }
}

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

class OpenChannelDNS {
  constructor(opt) {
    const o = Object.assign({
      Re: 2000, Lx: Math.PI, Lz: Math.PI / 2, Nx: 32, Ny: 48, Nz: 32,
      stretch: 1.8, cfl: 1.2, init: "turbulent", amp: 0.3, seed: 1, tavg: 50,
    }, opt || {});
    this.opt = o;
    this.Re = o.Re; this.nu = 1 / o.Re;
    this.Lx = o.Lx; this.Lz = o.Lz;
    const Nx = (this.Nx = o.Nx), Ny = (this.Ny = o.Ny), Nz = (this.Nz = o.Nz);
    this.dx = o.Lx / Nx; this.dz = o.Lz / Nz;
    this.cfl = o.cfl;
    const n = (this.n = Nx * Nz);

    // --- wall-normal grid (clustered at the wall y=0 only) ---
    const yf = (this.yf = new Float64Array(Ny + 1));
    for (let j = 0; j <= Ny; j++) {
      const eta = j / Ny;
      yf[j] = o.stretch > 1e-6 ? 1 - Math.tanh(o.stretch * (1 - eta)) / Math.tanh(o.stretch) : eta;
    }
    yf[0] = 0; yf[Ny] = 1;
    const yc = (this.yc = new Float64Array(Ny));
    const dyc = (this.dyc = new Float64Array(Ny)); // cell heights
    for (let j = 0; j < Ny; j++) { yc[j] = 0.5 * (yf[j] + yf[j + 1]); dyc[j] = yf[j + 1] - yf[j]; }
    const dyf = (this.dyf = new Float64Array(Ny + 1)); // centre-to-centre distances at faces
    for (let j = 1; j < Ny; j++) dyf[j] = yc[j] - yc[j - 1];
    dyf[0] = 2 * yc[0];            // mirror ghost below the wall
    dyf[Ny] = 2 * (1 - yc[Ny - 1]); // mirror ghost above the surface

    // periodic neighbour tables
    this.ip = new Int32Array(Nx); this.im = new Int32Array(Nx);
    for (let i = 0; i < Nx; i++) { this.ip[i] = (i + 1) % Nx; this.im[i] = (i - 1 + Nx) % Nx; }
    this.kp = new Int32Array(Nz); this.km = new Int32Array(Nz);
    for (let k = 0; k < Nz; k++) { this.kp[k] = (k + 1) % Nz; this.km[k] = (k - 1 + Nz) % Nz; }

    // --- fields: index = i + Nx*k + n*j ---
    const N3 = n * Ny, N3v = n * (Ny + 1);
    this.u = new Float64Array(N3); this.w = new Float64Array(N3); this.v = new Float64Array(N3v);
    this.Hu = new Float64Array(N3); this.Hw = new Float64Array(N3); this.Hv = new Float64Array(N3v);
    this.Hu1 = new Float64Array(N3); this.Hw1 = new Float64Array(N3); this.Hv1 = new Float64Array(N3v);
    this.du = new Float64Array(N3); this.dw = new Float64Array(N3); this.dv = new Float64Array(N3v);
    this.phi = new Float64Array(N3);
    this.pr = new Float64Array(N3); this.pi = new Float64Array(N3); this.cp = new Float64Array(N3);

    // Poisson operators
    this.fx = new FFT1(Nx); this.fz = new FFT1(Nz);
    this.K2 = new Float64Array(n);
    for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
      const kx = (2 / this.dx) * Math.sin((Math.PI * i) / Nx);
      const kz = (2 / this.dz) * Math.sin((Math.PI * k) / Nz);
      this.K2[i + Nx * k] = kx * kx + kz * kz;
    }
    this.pa = new Float64Array(Ny); this.pc = new Float64Array(Ny);
    for (let j = 0; j < Ny; j++) {
      this.pa[j] = j > 0 ? 1 / (dyc[j] * dyf[j]) : 0;
      this.pc[j] = j < Ny - 1 ? 1 / (dyc[j] * dyf[j + 1]) : 0;
    }
    // implicit-diffusion scratch
    this.tcp = new Float64Array(Ny + 1); this.tden = new Float64Array(Ny + 1);

    // RK3 (Wray) coefficients
    this.G = [8 / 15, 5 / 12, 3 / 4];
    this.Z = [0, -17 / 60, -5 / 12];
    this.A = [8 / 15, 2 / 15, 1 / 3];

    this.t = 0; this.nstep = 0; this.dt = 0;
    this.dtVisc = 2.0 / (4 * this.nu * (1 / (this.dx * this.dx) + 1 / (this.dz * this.dz)));
    this.dpdx = 0;
    this.history = [];
    this.resetStats();
    this.initialize();
  }

  // ---------- initial condition ----------
  initialize() {
    const { Nx, Ny, Nz, n, u, v, w, yc, yf, dx, dz, Lx, Lz } = this;
    const o = this.opt, rnd = mulberry32(o.seed | 0 || 1), A = o.amp;
    const U = (y) => (o.init === "laminar" ? 1.5 * (2 * y - y * y) : (8 / 7) * Math.pow(Math.max(y, 0), 1 / 7));
    // large-scale modes: (streamwise index, spanwise index)
    const modes = [[0, 1], [0, 2], [1, 1], [1, 2], [2, 1], [1, 3], [2, 2], [1, 0]];
    const comp = () => modes.map(([mx, mz]) => ({
      a: (2 * Math.PI * mx) / Lx, b: (2 * Math.PI * mz) / Lz,
      p1: 2 * Math.PI * rnd(), p2: 2 * Math.PI * rnd(), c: (0.5 + rnd()) / Math.sqrt(modes.length),
    }));
    const mu = comp(), mv = comp(), mw = comp();
    const field = (ms, x, z) => {
      let s = 0;
      for (const m of ms) s += m.c * Math.cos(m.a * x + m.p1) * Math.cos(m.b * z + m.p2);
      return s;
    };
    for (let j = 0; j < Ny; j++) {
      const env = Math.sin(Math.PI * yc[j]);
      for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
        const I = i + Nx * k + n * j;
        const xc = (i + 0.5) * dx, zc = (k + 0.5) * dz;
        u[I] = U(yc[j]) + A * env * (field(mu, i * dx, zc) + 0.1 * (rnd() - 0.5));
        w[I] = A * env * (field(mw, xc, k * dz) + 0.1 * (rnd() - 0.5));
      }
    }
    for (let j = 1; j < Ny; j++) {
      const env = Math.sin(Math.PI * yf[j]);
      for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
        v[i + Nx * k + n * j] = A * env * (field(mv, (i + 0.5) * dx, (k + 0.5) * dz) + 0.1 * (rnd() - 0.5));
      }
    }
    // remove any plane-mean spanwise velocity introduced by the random noise
    for (let j = 0; j < Ny; j++) {
      const mw = this.planeMean(w, j), o = n * j;
      for (let m = 0; m < n; m++) w[o + m] -= mw;
    }
    this.project(1);
    this.fixFlux(1);
    this.computeDt(true);
    this.record();
  }

  // ---------- right-hand sides ----------
  // H = -convection + nu*(d2/dx2 + d2/dz2);  L = nu*d2/dy2 (stored into du/dv/dw)
  computeRHS() {
    const { Nx, Ny, Nz, n, u, v, w, Hu, Hv, Hw, du, dv, dw, ip, im, kp, km, dyc, dyf } = this;
    const nu = this.nu, idx = 1 / this.dx, idz = 1 / this.dz;
    const nx2 = nu * idx * idx, nz2 = nu * idz * idz;
    for (let j = 0; j < Ny; j++) {
      const oj = n * j, ojp = j < Ny - 1 ? oj + n : oj, ojm = j > 0 ? oj - n : oj;
      const vT = n * (j + 1), vB = n * j; // v levels at top/bottom face of u/w cell j
      const ihy = 1 / dyc[j];
      const aLy = nu / (dyc[j] * dyf[j]), aUy = nu / (dyc[j] * dyf[j + 1]);
      const top = j < Ny - 1, bot = j > 0;
      for (let k = 0; k < Nz; k++) {
        const kP = kp[k], kM = km[k];
        const rk = Nx * k, rkP = Nx * kP, rkM = Nx * kM;
        for (let i = 0; i < Nx; i++) {
          const iP = ip[i], iM = im[i];
          const I = i + rk + oj;
          // ----- u at (x_i face, yc_j, zc_k) -----
          const u0 = u[I], uE = u[iP + rk + oj], uW = u[iM + rk + oj];
          const uN = u[i + rkP + oj], uS = u[i + rkM + oj];
          let ucE = 0.5 * (u0 + uE), ucW = 0.5 * (uW + u0);
          let cu = (ucE * ucE - ucW * ucW) * idx;
          cu += (0.5 * (u0 + uN) * 0.5 * (w[iM + rkP + oj] + w[i + rkP + oj])
               - 0.5 * (u0 + uS) * 0.5 * (w[iM + rk + oj] + w[i + rk + oj])) * idz;
          const uT = u[i + rk + ojp], uB = u[i + rk + ojm];
          const fT = top ? 0.5 * (u0 + uT) * 0.5 * (v[iM + rk + vT] + v[i + rk + vT]) : 0;
          const fB = bot ? 0.5 * (u0 + uB) * 0.5 * (v[iM + rk + vB] + v[i + rk + vB]) : 0;
          cu += (fT - fB) * ihy;
          Hu[I] = -cu + nx2 * (uE - 2 * u0 + uW) + nz2 * (uN - 2 * u0 + uS);
          const ugT = top ? uT : u0, ugB = bot ? uB : -u0;
          du[I] = aUy * (ugT - u0) - aLy * (u0 - ugB);

          // ----- w at (xc_i, yc_j, z_k face) -----
          const w0 = w[I], wN = w[i + rkP + oj], wS = w[i + rkM + oj];
          const wE = w[iP + rk + oj], wW = w[iM + rk + oj];
          const wcN = 0.5 * (w0 + wN), wcS = 0.5 * (wS + w0);
          let cw = (wcN * wcN - wcS * wcS) * idz;
          cw += (0.5 * (u[iP + rkM + oj] + u[iP + rk + oj]) * 0.5 * (w0 + wE)
               - 0.5 * (u[i + rkM + oj] + u0) * 0.5 * (wW + w0)) * idx;
          const wT = w[i + rk + ojp], wB = w[i + rk + ojm];
          const gT = top ? 0.5 * (w0 + wT) * 0.5 * (v[i + rkM + vT] + v[i + rk + vT]) : 0;
          const gB = bot ? 0.5 * (w0 + wB) * 0.5 * (v[i + rkM + vB] + v[i + rk + vB]) : 0;
          cw += (gT - gB) * ihy;
          Hw[I] = -cw + nx2 * (wE - 2 * w0 + wW) + nz2 * (wN - 2 * w0 + wS);
          const wgT = top ? wT : w0, wgB = bot ? wB : -w0;
          dw[I] = aUy * (wgT - w0) - aLy * (w0 - wgB);
        }
      }
    }
    // ----- v at (xc_i, y_j face, zc_k), j = 1..Ny-1 -----
    for (let j = 1; j < Ny; j++) {
      const oj = n * j, ojm = oj - n; // u/w cells j (above) and j-1 (below)
      const ihy = 1 / dyf[j];
      const bL = nu / (dyf[j] * dyc[j - 1]), bU = nu / (dyf[j] * dyc[j]);
      for (let k = 0; k < Nz; k++) {
        const kP = kp[k], kM = km[k];
        const rk = Nx * k, rkP = Nx * kP, rkM = Nx * kM;
        for (let i = 0; i < Nx; i++) {
          const iP = ip[i], iM = im[i];
          const I = i + rk + oj;
          const v0 = v[I], vE = v[iP + rk + oj], vW = v[iM + rk + oj];
          const vN = v[i + rkP + oj], vS = v[i + rkM + oj];
          const vT = v[I + n], vB = v[I - n];
          let c = (0.5 * (u[iP + rk + ojm] + u[iP + rk + oj]) * 0.5 * (v0 + vE)
                 - 0.5 * (u[i + rk + ojm] + u[i + rk + oj]) * 0.5 * (vW + v0)) * idx;
          const vcT = 0.5 * (v0 + vT), vcB = 0.5 * (vB + v0);
          c += (vcT * vcT - vcB * vcB) * ihy;
          c += (0.5 * (w[i + rkP + ojm] + w[i + rkP + oj]) * 0.5 * (v0 + vN)
              - 0.5 * (w[i + rk + ojm] + w[i + rk + oj]) * 0.5 * (vS + v0)) * idz;
          Hv[I] = -c + nx2 * (vE - 2 * v0 + vW) + nz2 * (vN - 2 * v0 + vS);
          dv[I] = bU * (vT - v0) - bL * (v0 - vB);
        }
      }
    }
  }

  // ---------- implicit wall-normal diffusion: (I - beta*L) x = rhs ----------
  solveImplicitUW(arr, beta) {
    const { Ny, n, dyc, dyf, tcp, tden } = this, nu = this.nu;
    // build and factorise (coefficients are the same for every (i,k))
    let prevC = 0;
    for (let j = 0; j < Ny; j++) {
      const aL = nu / (dyc[j] * dyf[j]), aU = nu / (dyc[j] * dyf[j + 1]);
      let lo = 0, di, up = 0;
      if (Ny === 1) di = 1 + beta * 2 * aL;
      else if (j === 0) { di = 1 + beta * (2 * aL + aU); up = -beta * aU; }
      else if (j === Ny - 1) { di = 1 + beta * aL; lo = -beta * aL; }
      else { di = 1 + beta * (aL + aU); lo = -beta * aL; up = -beta * aU; }
      const den = di - lo * prevC;
      tden[j] = den; tcp[j] = up / den; prevC = tcp[j];
      this._lo = this._lo || new Float64Array(Ny + 1); this._lo[j] = lo;
    }
    const lo = this._lo;
    for (let j = 0; j < Ny; j++) {
      const o = n * j, om = o - n, l = lo[j], id = 1 / tden[j];
      if (j === 0) for (let m = 0; m < n; m++) arr[m] *= id;
      else for (let m = 0; m < n; m++) arr[o + m] = (arr[o + m] - l * arr[om + m]) * id;
    }
    for (let j = Ny - 2; j >= 0; j--) {
      const o = n * j, op = o + n, c = tcp[j];
      for (let m = 0; m < n; m++) arr[o + m] -= c * arr[op + m];
    }
  }

  solveImplicitV(arr, beta) {
    const { Ny, n, dyc, dyf, tcp, tden } = this, nu = this.nu;
    if (Ny < 2) return;
    const lo = this._lov || (this._lov = new Float64Array(Ny + 1));
    let prevC = 0;
    for (let j = 1; j < Ny; j++) {
      const bL = nu / (dyf[j] * dyc[j - 1]), bU = nu / (dyf[j] * dyc[j]);
      const di = 1 + beta * (bL + bU);
      const l = j > 1 ? -beta * bL : 0, up = j < Ny - 1 ? -beta * bU : 0;
      const den = di - l * prevC;
      tden[j] = den; tcp[j] = up / den; prevC = tcp[j]; lo[j] = l;
    }
    for (let j = 1; j < Ny; j++) {
      const o = n * j, om = o - n, l = lo[j], id = 1 / tden[j];
      if (j === 1) for (let m = 0; m < n; m++) arr[o + m] *= id;
      else for (let m = 0; m < n; m++) arr[o + m] = (arr[o + m] - l * arr[om + m]) * id;
    }
    for (let j = Ny - 2; j >= 1; j--) {
      const o = n * j, op = o + n, c = tcp[j];
      for (let m = 0; m < n; m++) arr[o + m] -= c * arr[op + m];
    }
  }

  // ---------- pressure projection ----------
  // Transform along one axis for `count` lines at once (inner loop over lines).
  fftAxis(re, im, off, sign, f, stride, count, cstride) {
    const N = f.N;
    if (N === 1) return;
    if (!f.pow2) {
      for (let c = 0; c < count; c++) {
        const o = off + c * cstride;
        for (let a = 0; a < N; a++) { f.re[a] = re[o + a * stride]; f.im[a] = im[o + a * stride]; }
        f.run(sign);
        for (let a = 0; a < N; a++) { re[o + a * stride] = f.re[a]; im[o + a * stride] = f.im[a]; }
      }
      return;
    }
    const rev = f.rev, cs = f.cs, sn = f.sn;
    for (let a = 0; a < N; a++) {
      const r = rev[a];
      if (r > a) {
        const oa = off + a * stride, ob = off + r * stride;
        for (let c = 0; c < count; c++) {
          const p = oa + c * cstride, q = ob + c * cstride;
          let t = re[p]; re[p] = re[q]; re[q] = t;
          t = im[p]; im[p] = im[q]; im[q] = t;
        }
      }
    }
    for (let size = 2; size <= N; size <<= 1) {
      const half = size >> 1, step = N / size;
      for (let k = 0; k < half; k++) {
        const wr = cs[k * step], wi = sign * sn[k * step];
        for (let start = 0; start < N; start += size) {
          const oa = off + (start + k) * stride, ob = oa + half * stride;
          for (let c = 0; c < count; c++) {
            const p = oa + c * cstride, q = ob + c * cstride;
            const qr = re[q], qi = im[q];
            const xr = qr * wr - qi * wi, xi = qr * wi + qi * wr;
            const pr_ = re[p], pi_ = im[p];
            re[q] = pr_ - xr; im[q] = pi_ - xi; re[p] = pr_ + xr; im[p] = pi_ + xi;
          }
        }
      }
    }
  }
  fft2(re, im, off, sign) {
    const { Nx, Nz } = this;
    this.fftAxis(re, im, off, sign, this.fx, 1, Nz, Nx);  // along x
    this.fftAxis(re, im, off, sign, this.fz, Nx, Nx, 1);  // along z
  }

  // Precompute the y-tridiagonal factorisation for every (kx,kz) (data independent).
  factorPoisson() {
    const { Ny, n, K2, pa, pc } = this;
    const cp = this.cp, inv = (this.pinv = new Float64Array(n * Ny));
    for (let j = 0; j < Ny; j++) {
      const o = n * j, om = o - n, a = pa[j], c = pc[j], b0 = -(pa[j] + pc[j]);
      for (let m = 0; m < n; m++) {
        if (m === 0 && j === 0) { cp[o] = 0; inv[o] = 0; continue; } // phi(k=0, j=0) = 0
        const den = b0 - K2[m] - (j > 0 ? a * cp[om + m] : 0);
        inv[o + m] = 1 / den; cp[o + m] = c / den;
      }
    }
    // mirror index (-kx,-kz) for packing two real planes into one complex FFT
    const { Nx, Nz } = this;
    this.mir = new Int32Array(n);
    for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++)
      this.mir[i + Nx * k] = ((Nx - i) % Nx) + Nx * ((Nz - k) % Nz);
    this.tr2 = new Float64Array(n); this.ti2 = new Float64Array(n);
  }

  // Make the velocity field discretely divergence-free; adt = alpha*dt.
  project(adt) {
    const { Nx, Ny, Nz, n, u, v, w, ip, kp, dyc, dyf, pr, pi, cp, pa, phi } = this;
    if (!this.pinv) this.factorPoisson();
    const inv = this.pinv, mir = this.mir, tr = this.tr2, ti = this.ti2;
    const idx = 1 / this.dx, idz = 1 / this.dz, s = 1 / adt;
    // divergence -> pr
    for (let j = 0; j < Ny; j++) {
      const oj = n * j, ihy = 1 / dyc[j];
      for (let k = 0; k < Nz; k++) {
        const rk = Nx * k, rkP = Nx * kp[k];
        for (let i = 0; i < Nx; i++) {
          const I = i + rk + oj;
          const d = (u[ip[i] + rk + oj] - u[I]) * idx + (v[I + n] - v[I]) * ihy + (w[i + rkP + oj] - w[I]) * idz;
          pr[I] = d * s; pi[I] = 0;
        }
      }
    }
    // forward FFTs, two real planes per complex transform
    for (let j = 0; j < Ny; j += 2) {
      const o = n * j;
      if (j + 1 >= Ny) { this.fft2(pr, pi, o, -1); break; }
      const o2 = o + n;
      for (let m = 0; m < n; m++) pi[o + m] = pr[o2 + m];
      this.fft2(pr, pi, o, -1);
      for (let m = 0; m < n; m++) { tr[m] = pr[o + m]; ti[m] = pi[o + m]; }
      for (let m = 0; m < n; m++) {
        const zr = tr[m], zi = ti[m], mr = tr[mir[m]], mi = -ti[mir[m]]; // conj Z(-k)
        pr[o + m] = 0.5 * (zr + mr); pi[o + m] = 0.5 * (zi + mi);
        pr[o2 + m] = 0.5 * (zi - mi); pi[o2 + m] = -0.5 * (zr - mr);
      }
    }
    // tridiagonal solve in y for every (kx,kz)
    for (let j = 0; j < Ny; j++) {
      const o = n * j, om = o - n, a = pa[j];
      if (j === 0) for (let m = 0; m < n; m++) { pr[m] *= inv[m]; pi[m] *= inv[m]; }
      else for (let m = 0; m < n; m++) {
        const id = inv[o + m];
        pr[o + m] = (pr[o + m] - a * pr[om + m]) * id; pi[o + m] = (pi[o + m] - a * pi[om + m]) * id;
      }
    }
    for (let j = Ny - 2; j >= 0; j--) {
      const o = n * j, op = o + n;
      for (let m = 0; m < n; m++) { const c = cp[o + m]; pr[o + m] -= c * pr[op + m]; pi[o + m] -= c * pi[op + m]; }
    }
    // inverse FFTs, two planes per transform
    const sc = 1 / n;
    for (let j = 0; j < Ny; j += 2) {
      const o = n * j;
      if (j + 1 >= Ny) {
        this.fft2(pr, pi, o, +1);
        for (let m = 0; m < n; m++) phi[o + m] = pr[o + m] * sc;
        break;
      }
      const o2 = o + n;
      for (let m = 0; m < n; m++) {
        const ar = pr[o + m], ai = pi[o + m], br = pr[o2 + m], bi = pi[o2 + m];
        pr[o + m] = ar - bi; pi[o + m] = ai + br;
      }
      this.fft2(pr, pi, o, +1);
      for (let m = 0; m < n; m++) { phi[o + m] = pr[o + m] * sc; phi[o2 + m] = pi[o + m] * sc; }
    }
    // correct velocities
    for (let j = 0; j < Ny; j++) {
      const oj = n * j;
      for (let k = 0; k < Nz; k++) {
        const rk = Nx * k, rkM = Nx * this.km[k];
        for (let i = 0; i < Nx; i++) {
          const I = i + rk + oj;
          u[I] -= adt * (phi[I] - phi[this.im[i] + rk + oj]) * idx;
          w[I] -= adt * (phi[I] - phi[i + rkM + oj]) * idz;
        }
      }
    }
    for (let j = 1; j < Ny; j++) {
      const oj = n * j, ihy = 1 / dyf[j];
      for (let m = 0; m < n; m++) v[oj + m] -= adt * (phi[oj + m] - phi[oj - n + m]) * ihy;
    }
  }

  // Bulk (volume-averaged) velocity of a cell-height-staggered component (u or w).
  bulk(arr) {
    const { Ny, n, dyc } = this;
    const a = arr || this.u;
    let ub = 0;
    for (let j = 0; j < Ny; j++) {
      let s = 0; const o = n * j;
      for (let m = 0; m < n; m++) s += a[o + m];
      ub += (s / n) * dyc[j];
    }
    return ub;
  }

  // Constant streamwise mass flux U_b = 1, applied as a uniform mean pressure gradient in x.
  // There is no mean spanwise pressure gradient: the spanwise flow is left free.
  // Returns -dP/dx of this substage.
  fixFlux(adt) {
    const cu = 1 - this.bulk(this.u);
    const u = this.u;
    for (let I = 0; I < u.length; I++) u[I] += cu;
    return cu / adt;
  }

  // Plane-averaged spanwise velocity: bulk value and largest |<w>(y)|.
  meanSpanwise() {
    let bulk = 0, peak = 0;
    for (let j = 0; j < this.Ny; j++) {
      const m = this.planeMean(this.w, j);
      bulk += m * this.dyc[j]; peak = Math.max(peak, Math.abs(m));
    }
    return { bulk, peak };
  }

  computeDt(init) {
    const { Nx, Ny, Nz, n, u, v, w, dyc } = this;
    const idx = 1 / this.dx, idz = 1 / this.dz;
    let mx = 1e-12;
    for (let j = 0; j < Ny; j++) {
      const oj = n * j, ihy = 1 / dyc[j];
      for (let m = 0; m < n; m++) {
        const I = oj + m;
        const vv = Math.max(Math.abs(v[I]), Math.abs(v[I + n]));
        const c = Math.abs(u[I]) * idx + vv * ihy + Math.abs(w[I]) * idz;
        if (c > mx) mx = c;
      }
    }
    this.cmax = mx;
    let dt = Math.min(this.cfl / mx, this.dtVisc);
    if (!init && this.dt > 0) dt = Math.min(dt, 1.1 * this.dt);
    if (init) dt *= 0.5;
    this.dt = dt;
    return dt;
  }

  step() {
    const dt = this.computeDt(false);
    this._idx = 1 / this.dx; this._idz = 1 / this.dz;
    let dp = 0;
    for (let s = 0; s < 3; s++) {
      const g = this.G[s], z = this.Z[s], a = this.A[s];
      this.computeRHS();
      const { u, v, w, Hu, Hv, Hw, Hu1, Hv1, Hw1, du, dv, dw } = this;
      for (let I = 0; I < u.length; I++) {
        du[I] = dt * (g * Hu[I] + z * Hu1[I] + a * du[I]);
        dw[I] = dt * (g * Hw[I] + z * Hw1[I] + a * dw[I]);
      }
      const n = this.n, nv = n * this.Ny;
      for (let I = n; I < nv; I++) dv[I] = dt * (g * Hv[I] + z * Hv1[I] + a * dv[I]);
      const beta = 0.5 * a * dt;
      this.solveImplicitUW(du, beta);
      this.solveImplicitUW(dw, beta);
      this.solveImplicitV(dv, beta);
      for (let I = 0; I < u.length; I++) { u[I] += du[I]; w[I] += dw[I]; }
      for (let I = n; I < nv; I++) v[I] += dv[I];
      // swap H buffers
      this.Hu = Hu1; this.Hu1 = Hu; this.Hv = Hv1; this.Hv1 = Hv; this.Hw = Hw1; this.Hw1 = Hw;
      this.project(a * dt);
      dp += a * dt * this.fixFlux(a * dt);
    }
    this.dpdx = dp / dt;
    this.t += dt;
    this.nstep++;
    this.accumulate(dt);
    if (!this.avgOn && this.t >= this.opt.tavg) { this.resetStats(); this.avgOn = true; }
    if (this.nstep % 5 === 0) this.record();
  }

  // ---------- diagnostics ----------
  planeMean(arr, j) {
    const n = this.n, o = n * j; let s = 0;
    for (let m = 0; m < n; m++) s += arr[o + m];
    return s / n;
  }
  wallShear() { return (this.nu * this.planeMean(this.u, 0)) / this.yc[0]; }
  reTau() { return Math.sqrt(Math.max(this.wallShear(), 0)) / this.nu; }
  maxDiv() {
    const { Nx, Ny, Nz, n, u, v, w, ip, kp, dyc } = this;
    let mx = 0;
    for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
      const I = i + Nx * k + n * j;
      const d = (u[ip[i] + Nx * k + n * j] - u[I]) / this.dx + (v[I + n] - v[I]) / dyc[j] + (w[i + Nx * kp[k] + n * j] - w[I]) / this.dz;
      mx = Math.max(mx, Math.abs(d));
    }
    return mx;
  }
  // turbulent kinetic energy of fluctuations (volume average)
  tke() {
    const { Ny, n, u, v, w, dyc } = this;
    let e = 0;
    for (let j = 0; j < Ny; j++) {
      const o = n * j, mu = this.planeMean(u, j), mw = this.planeMean(w, j);
      let s = 0;
      for (let m = 0; m < n; m++) {
        const a = u[o + m] - mu, c = w[o + m] - mw, b = 0.5 * (v[o + m] + v[o + n + m]);
        s += a * a + b * b + c * c;
      }
      e += 0.5 * (s / n) * dyc[j];
    }
    return e;
  }
  record() {
    this.history.push({ t: this.t, reTau: this.reTau(), tke: this.tke() });
    if (this.history.length > 4000) {
      // thin out: keep every other sample
      this.history = this.history.filter((_, i) => i % 2 === 0);
    }
  }

  resetStats() {
    const Ny = this.Ny;
    this.st = {
      T: 0, t0: this.t || 0, U: new Float64Array(Ny), W: new Float64Array(Ny), uu: new Float64Array(Ny), vv: new Float64Array(Ny),
      ww: new Float64Array(Ny), uv: new Float64Array(Ny), tw: 0,
    };
  }
  accumulate(dt) {
    const { Ny, n, u, v, w, st } = this;
    for (let j = 0; j < Ny; j++) {
      const o = n * j;
      let su = 0, sv = 0, sw = 0, suu = 0, svv = 0, sww = 0, suv = 0;
      for (let m = 0; m < n; m++) {
        const a = u[o + m], b = 0.5 * (v[o + m] + v[o + n + m]), c = w[o + m];
        su += a; sv += b; sw += c; suu += a * a; svv += b * b; sww += c * c; suv += a * b;
      }
      su /= n; sv /= n; sw /= n;
      st.U[j] += dt * su;
      st.W[j] += dt * sw;
      st.uu[j] += dt * (suu / n - su * su);
      st.vv[j] += dt * (svv / n - sv * sv);
      st.ww[j] += dt * (sww / n - sw * sw);
      st.uv[j] += dt * (suv / n - su * sv);
    }
    st.tw += dt * this.wallShear();
    st.T += dt;
  }
  stats() {
    const { st, Ny } = this;
    if (st.T <= 0) return null;
    const T = st.T, tw = st.tw / T, ut = Math.sqrt(Math.max(tw, 1e-30));
    const out = { T, t0: st.t0, utau: ut, reTau: ut / this.nu, yp: [], Up: [], Wp: [], urms: [], vrms: [], wrms: [], uv: [] };
    for (let j = 0; j < Ny; j++) {
      out.yp.push((this.yc[j] * ut) / this.nu);
      out.Up.push(st.U[j] / T / ut);
      out.Wp.push(st.W[j] / T / ut);
      out.urms.push(Math.sqrt(Math.max(st.uu[j] / T, 0)) / ut);
      out.vrms.push(Math.sqrt(Math.max(st.vv[j] / T, 0)) / ut);
      out.wrms.push(Math.sqrt(Math.max(st.ww[j] / T, 0)) / ut);
      out.uv.push(-st.uv[j] / T / (ut * ut));
    }
    return out;
  }

  startAveragingNow() { this.resetStats(); this.avgOn = true; }

  // Averaged statistics once averaging has started, otherwise instantaneous plane averages.
  profiles() {
    const st = this.avgOn ? this.stats() : null;
    if (st) return Object.assign(st, { averaged: true });
    const ut = Math.sqrt(Math.max(this.wallShear(), 1e-12));
    const out = { averaged: false, T: 0, t0: 0, reTau: ut / this.nu, yp: [], Up: [], Wp: [], urms: [], vrms: [], wrms: [], uv: [] };
    const { Ny, n, u, v, w } = this;
    for (let j = 0; j < Ny; j++) {
      const o = n * j; let su = 0, sv = 0, sw = 0, suu = 0, svv = 0, sww = 0, suv = 0;
      for (let m = 0; m < n; m++) {
        const a = u[o + m], b = 0.5 * (v[o + m] + v[o + n + m]), c = w[o + m];
        su += a; sv += b; sw += c; suu += a * a; svv += b * b; sww += c * c; suv += a * b;
      }
      su /= n; sv /= n; sw /= n;
      out.yp.push((this.yc[j] * ut) / this.nu); out.Up.push(su / ut); out.Wp.push(sw / ut);
      out.urms.push(Math.sqrt(Math.max(suu / n - su * su, 0)) / ut);
      out.vrms.push(Math.sqrt(Math.max(svv / n - sv * sv, 0)) / ut);
      out.wrms.push(Math.sqrt(Math.max(sww / n - sw * sw, 0)) / ut);
      out.uv.push(-(suv / n - su * sv) / (ut * ut));
    }
    return out;
  }

  // ---------- display fields (codes: 0 u, 1 u', 2 v, 3 w, 4 omega_z') ----------
  sampler(field) {
    const { Nx, Ny, n, u, v, w, ip, im, kp, yc } = this;
    const mean = new Float64Array(Ny);
    if (field === 1 || field === 4) for (let j = 0; j < Ny; j++) mean[j] = this.planeMean(u, j);
    const uc = (i, j, k) => 0.5 * (u[i + Nx * k + n * j] + u[ip[i] + Nx * k + n * j]);
    const wc = (i, j, k) => 0.5 * (w[i + Nx * k + n * j] + w[i + Nx * kp[k] + n * j]);
    const vc = (i, j, k) => 0.5 * (v[i + Nx * k + n * j] + v[i + Nx * k + n * (j + 1)]);
    switch (field) {
      case 0: return uc;
      case 1: return (i, j, k) => uc(i, j, k) - mean[j];
      case 2: return vc;
      case 3: return wc;
      default: // fluctuating spanwise vorticity omega_z' = dv'/dx - du'/dy (u = 0 at the wall, du/dy = 0 at the surface)
        return (i, j, k) => {
          const top = j < Ny - 1, bot = j > 0;
          const uT = top ? uc(i, j + 1, k) - mean[j + 1] : uc(i, j, k) - mean[j], yT = top ? yc[j + 1] : 2 - yc[j];
          const uB = bot ? uc(i, j - 1, k) - mean[j - 1] : -(uc(i, j, k) - mean[j]), yB = bot ? yc[j - 1] : -yc[j];
          return (vc(ip[i], j, k) - vc(im[i], j, k)) / (2 * this.dx) - (uT - uB) / (yT - yB);
        };
    }
  }
  // plan: [k*Nx + i] at y index j; side: [j*Nx + i] at z index k; end: [j*Nz + k] at x index i
  slices(field, si, sj, sk) {
    const { Nx, Ny, Nz } = this, f = this.sampler(field);
    const plan = new Float32Array(Nx * Nz), side = new Float32Array(Nx * Ny), end = new Float32Array(Nz * Ny);
    for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) plan[k * Nx + i] = f(i, sj, k);
    for (let j = 0; j < Ny; j++) for (let i = 0; i < Nx; i++) side[j * Nx + i] = f(i, j, sk);
    for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k++) end[j * Nz + k] = f(si, j, k);
    return { plan, side, end };
  }
  // Colour range: [0, max] for u, otherwise +-3 rms over y < 0.6
  fieldRange(field) {
    const { Nx, Ny, Nz } = this, f = this.sampler(field);
    if (field === 0) {
      let mx = 0;
      for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k += 2) for (let i = 0; i < Nx; i += 2) mx = Math.max(mx, f(i, j, k));
      return [0, Math.max(0.5, Math.ceil(mx * 10) / 10)];
    }
    let ss = 0, c = 0;
    for (let j = 0; j < Ny && this.yc[j] < 0.6; j++) for (let k = 0; k < Nz; k += 2) for (let i = 0; i < Nx; i += 2) { const a = f(i, j, k); ss += a * a; c++; }
    const lim = Number((3 * Math.sqrt(ss / Math.max(c, 1)) || 0.1).toPrecision(2));
    return [-lim, lim];
  }
}


if (typeof module !== "undefined") module.exports = { OpenChannelDNS, FFT1 };
