const { OpenChannelDNS, FFT1 } = require('../src/solver.js');
// FFT check for non-pow2
for (const N of [16, 24]) { const f = new FFT1(N); for (let i=0;i<N;i++){f.re[i]=Math.sin(3*2*Math.PI*i/N);f.im[i]=0;} f.run(-1); f.run(1); let e=0; for(let i=0;i<N;i++) e=Math.max(e,Math.abs(f.re[i]/N-Math.sin(3*2*Math.PI*i/N))); console.log('fft',N,e.toExponential(2)); }
const opts = JSON.parse(process.argv[2] || '{}');
const s = new OpenChannelDNS(opts);
console.log('init div', s.maxDiv().toExponential(2), 'Ub', s.bulk().toFixed(6), 'yplus1', (s.yc[0]*s.reTau()).toFixed(2));
const T = +(process.argv[3] || 60);
let t0 = Date.now(), last = 0;
while (s.t < T) {
  s.step();
  if (s.t - last > T/15) { last = s.t; console.log(`t=${s.t.toFixed(1)} dt=${s.dt.toFixed(4)} Re_tau=${s.reTau().toFixed(1)} tke=${s.tke().toFixed(5)} div=${s.maxDiv().toExponential(1)} Ub=${s.bulk().toFixed(6)} -dPdx=${s.dpdx.toFixed(5)}`); }
  if (!isFinite(s.u[0])) { console.log('NaN'); break; }
}
console.log('steps', s.nstep, 'ms/step', ((Date.now()-t0)/s.nstep).toFixed(2));
