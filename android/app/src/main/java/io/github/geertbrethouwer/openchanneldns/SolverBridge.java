package io.github.geertbrethouwer.openchanneldns;

import android.util.Base64;
import android.webkit.JavascriptInterface;

import org.json.JSONException;
import org.json.JSONObject;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;

/**
 * Exposed to the page as window.OCDNSNative. The page's NativeBackend calls these methods;
 * small results are returned as JSON, slice data as base64 of little-endian float32.
 */
final class SolverBridge {
    private volatile boolean wantRunning = false;
    private volatile boolean suspended = false;
    private volatile boolean created = false;

    @JavascriptInterface
    public void init(String json) throws JSONException {
        JSONObject p = new JSONObject(json);
        double[] a = {
            p.getDouble("Re"), p.getDouble("Lx"), p.getDouble("Lz"),
            p.getInt("Nx"), p.getInt("Ny"), p.getInt("Nz"),
            p.getDouble("stretch"), p.getDouble("cfl"), p.getDouble("amp"),
            p.optDouble("tavg", 50), "laminar".equals(p.optString("init")) ? 1 : 0, p.optInt("seed", 1),
            Math.max(1, Math.min(p.optInt("threads", 4), cores()))
        };
        NativeSolver.create(a);
        created = true;
        apply();
    }

    /** Number of CPU cores the page may offer as thread counts. */
    @JavascriptInterface
    public int cores() {
        return Runtime.getRuntime().availableProcessors();
    }

    @JavascriptInterface
    public void setRunning(boolean on) {
        wantRunning = on;
        apply();
    }

    /** Called by the activity: pause the solver while the app is in the background. */
    void setSuspended(boolean s) {
        suspended = s;
        apply();
    }

    void destroy() {
        if (created) NativeSolver.destroy();
        created = false;
    }

    private void apply() {
        if (created) NativeSolver.setRunning(wantRunning && !suspended);
    }

    @JavascriptInterface
    public String grid() {
        double[] g = NativeSolver.grid();
        int ny = (int) g[1];
        StringBuilder sb = new StringBuilder(64 + 24 * (2 * ny + 1));
        sb.append("{\"Nx\":").append((int) g[0]).append(",\"Ny\":").append(ny).append(",\"Nz\":").append((int) g[2])
          .append(",\"dx\":").append(num(g[3])).append(",\"dz\":").append(num(g[4]))
          .append(",\"threads\":").append((int) g[5]).append(",\"yc\":");
        array(sb, g, 6, ny);
        sb.append(",\"yf\":");
        array(sb, g, 6 + ny, ny + 1);
        return sb.append('}').toString();
    }

    @JavascriptInterface
    public String state() {
        double[] s = NativeSolver.state();
        return "{\"t\":" + num(s[0]) + ",\"dt\":" + num(s[1]) + ",\"nstep\":" + (long) s[2]
            + ",\"reTau\":" + num(s[3]) + ",\"blown\":" + (s[4] != 0) + "}";
    }

    @JavascriptInterface
    public String slices(int field, int i, int j, int k) {
        float[] a = NativeSolver.slices(field, i, j, k);
        ByteBuffer bb = ByteBuffer.allocate(4 * a.length).order(ByteOrder.LITTLE_ENDIAN);
        bb.asFloatBuffer().put(a);
        return Base64.encodeToString(bb.array(), Base64.NO_WRAP);
    }

    @JavascriptInterface
    public String range(int field) {
        double[] r = NativeSolver.range(field);
        return "[" + num(r[0]) + "," + num(r[1]) + "]";
    }

    @JavascriptInterface
    public String profiles() {
        double[] p = NativeSolver.profiles();
        int ny = (p.length - 4) / 7;
        StringBuilder sb = new StringBuilder(128 + 7 * 24 * ny);
        sb.append("{\"averaged\":").append(p[0] != 0).append(",\"T\":").append(num(p[1]))
          .append(",\"t0\":").append(num(p[2])).append(",\"reTau\":").append(num(p[3]));
        String[] names = {"yp", "Up", "Wp", "urms", "vrms", "wrms", "uv"};
        for (int q = 0; q < names.length; q++) {
            sb.append(",\"").append(names[q]).append("\":");
            array(sb, p, 4 + q * ny, ny);
        }
        return sb.append('}').toString();
    }

    @JavascriptInterface
    public String history() {
        double[] h = NativeSolver.history();
        int n = (int) h[0];
        StringBuilder sb = new StringBuilder(32 + 3 * 20 * n);
        sb.append("{\"t\":");
        array(sb, h, 1, n);
        sb.append(",\"reTau\":");
        array(sb, h, 1 + n, n);
        sb.append(",\"tke\":");
        array(sb, h, 1 + 2 * n, n);
        return sb.append('}').toString();
    }

    @JavascriptInterface
    public void startAveragingNow() {
        NativeSolver.startAveragingNow();
    }

    private static String num(double x) {
        return Double.isFinite(x) ? Double.toString(x) : "null";
    }

    private static void array(StringBuilder sb, double[] a, int off, int n) {
        sb.append('[');
        for (int q = 0; q < n; q++) {
            if (q > 0) sb.append(',');
            sb.append(num(a[off + q]));
        }
        sb.append(']');
    }
}
