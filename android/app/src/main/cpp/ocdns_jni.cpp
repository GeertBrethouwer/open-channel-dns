// JNI bridge between the Android app and the C++ solver (native/ocdns.hpp).
// The solver steps continuously on its own thread while running; the UI takes
// short snapshots (slices, profiles, history) under the same mutex.
#include <jni.h>

#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "ocdns.hpp"

namespace {

class Runner {
 public:
  explicit Runner(const ocdns::Params& p) : sim(p), worker([this] { loop(); }) {}
  ~Runner() {
    {
      std::lock_guard<std::mutex> lk(m);
      quit = true;
    }
    cv.notify_all();
    worker.join();
  }
  void setRunning(bool on) {
    {
      std::lock_guard<std::mutex> lk(m);
      running = on;
    }
    cv.notify_all();
  }
  // Snapshot access: tell the stepping thread to let go of the lock between steps.
  struct Access {
    Runner& r;
    std::unique_lock<std::mutex> lk;
    explicit Access(Runner& rr) : r(rr), lk((rr.waiting.fetch_add(1), rr.m)) { r.waiting.fetch_sub(1); }
  };

  std::mutex m;
  std::condition_variable cv;
  std::atomic<int> waiting{0};
  ocdns::Solver sim;
  bool running = false, quit = false, blown = false;

 private:
  std::thread worker;
  void loop() {
    for (;;) {
      while (waiting.load() > 0) std::this_thread::yield();
      std::unique_lock<std::mutex> lk(m);
      cv.wait(lk, [&] { return quit || (running && !blown); });
      if (quit) return;
      sim.step();
      if (!sim.finite()) {
        blown = true;
        running = false;
      }
    }
  }
};

std::unique_ptr<Runner> g_runner;
std::mutex g_life;  // guards create/destroy against concurrent calls

jdoubleArray toJava(JNIEnv* env, const std::vector<double>& v) {
  jdoubleArray a = env->NewDoubleArray(jsize(v.size()));
  env->SetDoubleArrayRegion(a, 0, jsize(v.size()), v.data());
  return a;
}

}  // namespace

#define JNI_FN(name) Java_io_github_geertbrethouwer_openchanneldns_NativeSolver_##name

extern "C" {

JNIEXPORT void JNICALL JNI_FN(create)(JNIEnv* env, jclass, jdoubleArray params) {
  std::vector<double> a(13, 0.0);
  env->GetDoubleArrayRegion(params, 0, std::min<jsize>(13, env->GetArrayLength(params)), a.data());
  ocdns::Params p;
  p.Re = a[0]; p.Lx = a[1]; p.Lz = a[2];
  p.Nx = int(a[3]); p.Ny = int(a[4]); p.Nz = int(a[5]);
  p.stretch = a[6]; p.cfl = a[7]; p.amp = a[8]; p.tavg = a[9];
  p.laminar = a[10] != 0; p.seed = uint32_t(a[11]);
  p.threads = int(a[12]);  // 0 = default
  std::lock_guard<std::mutex> lk(g_life);
  g_runner.reset();
  g_runner = std::make_unique<Runner>(p);
}

JNIEXPORT void JNICALL JNI_FN(destroy)(JNIEnv*, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  g_runner.reset();
}

JNIEXPORT void JNICALL JNI_FN(setRunning)(JNIEnv*, jclass, jboolean on) {
  std::lock_guard<std::mutex> lk(g_life);
  if (g_runner) g_runner->setRunning(on);
}

JNIEXPORT jdoubleArray JNICALL JNI_FN(state)(JNIEnv* env, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return toJava(env, {0, 0, 0, 0, 0});
  Runner::Access acc(*g_runner);
  auto& s = g_runner->sim;
  return toJava(env, {s.t, s.dt, double(s.nstep), s.reTau(), g_runner->blown ? 1.0 : 0.0});
}

JNIEXPORT jdoubleArray JNICALL JNI_FN(grid)(JNIEnv* env, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return toJava(env, {});
  Runner::Access acc(*g_runner);
  auto& s = g_runner->sim;
  std::vector<double> g = {double(s.Nx), double(s.Ny), double(s.Nz), s.dx, s.dz, double(s.threads())};
  g.insert(g.end(), s.yc.begin(), s.yc.end());
  g.insert(g.end(), s.yf.begin(), s.yf.end());
  return toJava(env, g);
}

JNIEXPORT jfloatArray JNICALL JNI_FN(slices)(JNIEnv* env, jclass, jint field, jint i, jint j, jint k) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return env->NewFloatArray(0);
  Runner::Access acc(*g_runner);
  auto& s = g_runner->sim;
  const int Nx = s.Nx, Ny = s.Ny, Nz = s.Nz;
  i = std::min(std::max(i, 0), Nx - 1);
  j = std::min(std::max(j, 0), Ny - 1);
  k = std::min(std::max(k, 0), Nz - 1);
  std::vector<float> buf(size_t(Nx) * Nz + size_t(Nx) * Ny + size_t(Nz) * Ny);
  s.slices(field, i, j, k, buf.data(), buf.data() + size_t(Nx) * Nz, buf.data() + size_t(Nx) * Nz + size_t(Nx) * Ny);
  jfloatArray a = env->NewFloatArray(jsize(buf.size()));
  env->SetFloatArrayRegion(a, 0, jsize(buf.size()), buf.data());
  return a;
}

JNIEXPORT jdoubleArray JNICALL JNI_FN(range)(JNIEnv* env, jclass, jint field) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return toJava(env, {0, 1});
  Runner::Access acc(*g_runner);
  double lo, hi;
  g_runner->sim.fieldRange(field, lo, hi);
  return toJava(env, {lo, hi});
}

JNIEXPORT jdoubleArray JNICALL JNI_FN(profiles)(JNIEnv* env, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return toJava(env, {0, 0, 0, 0});
  Runner::Access acc(*g_runner);
  const ocdns::Profiles p = g_runner->sim.profiles();
  std::vector<double> out = {p.averaged ? 1.0 : 0.0, p.T, p.t0, p.reTau};
  for (const auto* v : {&p.yp, &p.Up, &p.Wp, &p.urms, &p.vrms, &p.wrms, &p.uv}) out.insert(out.end(), v->begin(), v->end());
  return toJava(env, out);
}

JNIEXPORT jdoubleArray JNICALL JNI_FN(history)(JNIEnv* env, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return toJava(env, {0});
  Runner::Access acc(*g_runner);
  const auto& h = g_runner->sim.history;
  std::vector<double> out(1 + 3 * h.size());
  out[0] = double(h.size());
  for (size_t q = 0; q < h.size(); ++q) {
    out[1 + q] = h[q].t;
    out[1 + h.size() + q] = h[q].reTau;
    out[1 + 2 * h.size() + q] = h[q].tke;
  }
  return toJava(env, out);
}

JNIEXPORT void JNICALL JNI_FN(startAveragingNow)(JNIEnv*, jclass) {
  std::lock_guard<std::mutex> lk(g_life);
  if (!g_runner) return;
  Runner::Access acc(*g_runner);
  g_runner->sim.startAveragingNow();
}

}  // extern "C"
