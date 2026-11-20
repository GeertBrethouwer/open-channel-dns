"use strict";
/*
 * Manufactured-solution (MMS) test of src/solver.js.
 *
 * Exact solution (h = 1, free-slip surface at y = 1, no-slip wall at y = 0), three-dimensional,
 * time dependent and divergence free; a = 2 pi/Lx, b = 2 pi/Lz:
 *   x-y mode   u =  A1 G'(y) sin(ax),             v = -a A1 G(y) cos(ax)
 *   z-y mode   w =  A2 G'(y) sin(bz),             v = -b A2 G(y) cos(bz)
 *   x-z mode   u = 2b A3 P(y) cos(ax + 2bz),      w = -a A3 P(y) cos(ax + 2bz)
 *   mean       u = U(y) = 1.5 (2y - y^2) (bulk 1), w = Wm(t) P(y)
 * with G = y^2 (1 - y)(3 - 2y) (G = G' = 0 at y = 0, G = G'' = 0 at y = 1) and P = sin(pi y/2)
 * (P = 0 at y = 0, P' = 0 at y = 1). The pressure is zero; the body force
 *   F = du/dt + (u . grad) u - nu lap u
 * is added to the explicit right-hand side at every Runge-Kutta stage time (solver.forcing).
 *
 * Run with Node.js:    node test/mms.js [test]
 *      or JavaScriptCore: jsc src/solver.js test/mms.js -- [test]
 * test = space (grid refinement), time (dt refinement), all (default)
 */
const out = typeof print === "function" ? print : console.log;
const DNS = typeof OpenChannelDNS !== "undefined" ? OpenChannelDNS : require("../src/solver.js").OpenChannelDNS;
const argv = typeof process !== "undefined" && process.argv ? process.argv.slice(2)
  : typeof scriptArgs !== "undefined" ? scriptArgs : typeof arguments !== "undefined" ? Array.from(arguments) : [];

// ---------- exact solution: each component is a sum of T(t) Y(y) X(theta), theta = kx x + kz z ----------
function exactModel(Lx, Lz) {
  const a = (2 * Math.PI) / Lx, b = (2 * Math.PI) / Lz, h = Math.PI / 2;
  const e1 = 0.3, e2 = 0.25, e3 = 0.15, e4 = 0.2;
  // time factors [f, f']
  const A1 = (t) => [e1 * (1 + 0.5 * Math.sin(2 * t)), e1 * Math.cos(2 * t)];
  const A2 = (t) => [e2 * Math.cos(1.5 * t), -1.5 * e2 * Math.sin(1.5 * t)];
  const A3 = (t) => [e3 * (1 + 0.5 * Math.cos(3 * t)), -1.5 * e3 * Math.sin(3 * t)];
  const Wm = (t) => [e4 * Math.sin(t), e4 * Math.cos(t)];
  const one = () => [1, 0];
  const scale = (c, T) => (t) => { const r = T(t); return [c * r[0], c * r[1]]; };
  // wall-normal factors [Y, Y', Y'']
  const G = (y) => [3 * y * y - 5 * y ** 3 + 2 * y ** 4, 6 * y - 15 * y * y + 8 * y ** 3, 6 - 30 * y + 24 * y * y];
  const Gp = (y) => [6 * y - 15 * y * y + 8 * y ** 3, 6 - 30 * y + 24 * y * y, -30 + 48 * y];
  const P = (y) => [Math.sin(h * y), h * Math.cos(h * y), -h * h * Math.sin(h * y)];
  const U = (y) => [1.5 * (2 * y - y * y), 3 - 3 * y, -3];
  const term = (T, Y, trig, kx, kz) => ({ T, Y, sin: trig === "sin", kx, kz });
  return {
    u: [term(A1, Gp, "sin", a, 0), term(scale(2 * b, A3), P, "cos", a, 2 * b), term(one, U, "cos", 0, 0)],
    v: [term(scale(-a, A1), G, "cos", a, 0), term(scale(-b, A2), G, "cos", 0, b)],
    w: [term(A2, Gp, "sin", 0, b), term(scale(-a, A3), P, "cos", a, 2 * b), term(Wm, P, "cos", 0, 0)],
  };
}

// value, d/dt, d/dx, d/dy, d/dz, laplacian of one component at (x, y, z, t)
function evalComp(terms, x, y, z, t) {
  const r = [0, 0, 0, 0, 0, 0];
  for (const m of terms) {
    const T = m.T(t), Y = m.Y(y), th = m.kx * x + m.kz * z;
    const X = m.sin ? Math.sin(th) : Math.cos(th), dX = m.sin ? Math.cos(th) : -Math.sin(th);
    const k2 = m.kx * m.kx + m.kz * m.kz;
    r[0] += T[0] * Y[0] * X;
    r[1] += T[1] * Y[0] * X;
    r[2] += T[0] * Y[0] * m.kx * dX;
    r[3] += T[0] * Y[1] * X;
    r[4] += T[0] * Y[0] * m.kz * dX;
    r[5] += T[0] * (Y[2] - k2 * Y[0]) * X;
  }
  return r;
}

function force(E, nu, x, y, z, t) {
  const u = evalComp(E.u, x, y, z, t), v = evalComp(E.v, x, y, z, t), w = evalComp(E.w, x, y, z, t);
  const f = (c) => c[1] + u[0] * c[2] + v[0] * c[3] + w[0] * c[4] - nu * c[5];
  return [f(u), f(v), f(w)];
}

// staggered positions: u (i dx, yc, (k+1/2) dz), w ((i+1/2) dx, yc, k dz), v ((i+1/2) dx, yf, (k+1/2) dz)
function setExact(s, E, t) {
  const { Nx, Ny, Nz, n, dx, dz, yc, yf } = s;
  for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
    const I = i + Nx * k + n * j;
    s.u[I] = evalComp(E.u, i * dx, yc[j], (k + 0.5) * dz, t)[0];
    s.w[I] = evalComp(E.w, (i + 0.5) * dx, yc[j], k * dz, t)[0];
  }
  s.v.fill(0);
  for (let j = 1; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++)
    s.v[i + Nx * k + n * j] = evalComp(E.v, (i + 0.5) * dx, yf[j], (k + 0.5) * dz, t)[0];
}

function errors(s, E, t) {
  const { Nx, Ny, Nz, n, dx, dz, yc, yf } = s;
  let eu = 0, ev = 0, ew = 0;
  for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
    const I = i + Nx * k + n * j;
    eu = Math.max(eu, Math.abs(s.u[I] - evalComp(E.u, i * dx, yc[j], (k + 0.5) * dz, t)[0]));
    ew = Math.max(ew, Math.abs(s.w[I] - evalComp(E.w, (i + 0.5) * dx, yc[j], k * dz, t)[0]));
    if (j > 0) ev = Math.max(ev, Math.abs(s.v[I] - evalComp(E.v, (i + 0.5) * dx, yf[j], (k + 0.5) * dz, t)[0]));
  }
  return [eu, ev, ew];
}

// Run the manufactured solution to time T with fixed dt; returns the final solver.
function run(opt, dt, T) {
  const s = new DNS(Object.assign({ Re: 100, Lx: 4, Lz: 2, stretch: 1.5, amp: 0 }, opt));
  const E = exactModel(s.Lx, s.Lz);
  const { Nx, Ny, Nz, n, dx, dz, yc, yf } = s;
  // force at the staggered points, cached per stage time
  s.forcing = (t, sv) => {
    for (let j = 0; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++) {
      const I = i + Nx * k + n * j;
      sv.Hu[I] += force(E, sv.nu, i * dx, yc[j], (k + 0.5) * dz, t)[0];
      sv.Hw[I] += force(E, sv.nu, (i + 0.5) * dx, yc[j], k * dz, t)[2];
    }
    for (let j = 1; j < Ny; j++) for (let k = 0; k < Nz; k++) for (let i = 0; i < Nx; i++)
      sv.Hv[i + Nx * k + n * j] += force(E, sv.nu, (i + 0.5) * dx, yf[j], (k + 0.5) * dz, t)[1];
  };
  s.t = 0;
  setExact(s, E, 0);
  s.project(1); s.fixFlux(1);
  s.Hu1.fill(0); s.Hv1.fill(0); s.Hw1.fill(0); s.dt = 0;
  const nst = Math.round(T / dt);
  for (let q = 0; q < nst; q++) s.step(dt);
  s.err = errors(s, E, s.t);
  return s;
}

const fmt = (x) => x.toExponential(2);
const which = argv[0] || "all";

if (which === "space" || which === "all") {
  for (const rk of [3, 4]) {
    out(`-- grid refinement, rk = ${rk}, Re = 100, dt = 1e-3, T = 0.5 (max error u, v, w; order)`);
    let prev = null;
    for (const N of [8, 16, 32]) {
      const s = run({ rk, Nx: N, Ny: N, Nz: N }, 1e-3, 0.5);
      const e = s.err, m = Math.max(...e);
      out(`N = ${String(N).padStart(3)}   ${e.map(fmt).join("  ")}   div ${fmt(s.maxDiv())}` +
          (prev ? `   order ${(Math.log(prev / m) / Math.LN2).toFixed(2)}` : ""));
      prev = m;
    }
  }
}

if (which === "time" || which === "all") {
  // error against a reference with a much smaller dt on the same grid (removes the spatial error)
  const grid = { Nx: 16, Ny: 16, Nz: 16 }, T = 1, Re = 1000;
  for (const rk of [3, 4]) {
    const ref = run(Object.assign({ rk, Re }, grid), 1 / 1280, T);
    out(`-- time refinement, rk = ${rk}, Re = ${Re}, 16^3, T = ${T}: max |u - u_ref| (u, v, w; order)`);
    let prev = null;
    for (const dt of [0.04, 0.02, 0.01, 0.005]) {
      const s = run(Object.assign({ rk, Re }, grid), dt, T);
      let m = 0;
      const e = ["u", "v", "w"].map((c) => { let x = 0; for (let I = 0; I < s[c].length; I++) x = Math.max(x, Math.abs(s[c][I] - ref[c][I])); m = Math.max(m, x); return x; });
      out(`dt = ${dt.toFixed(4)}   ${e.map(fmt).join("  ")}` + (prev ? `   order ${(Math.log(prev / m) / Math.LN2).toFixed(2)}` : ""));
      prev = m;
    }
  }
}

if (typeof module !== "undefined") module.exports = { exactModel, evalComp, force, run };
