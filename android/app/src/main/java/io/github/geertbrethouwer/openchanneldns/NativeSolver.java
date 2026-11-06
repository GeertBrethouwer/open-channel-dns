package io.github.geertbrethouwer.openchanneldns;

/** JNI entry points of the C++ solver (native/ocdns.hpp, single precision). One solver at a time. */
final class NativeSolver {
    static {
        System.loadLibrary("ocdns");
    }

    private NativeSolver() {}

    /** params: Re, Lx, Lz, Nx, Ny, Nz, stretch, cfl, amp, tavg, laminar (0/1), seed */
    static native void create(double[] params);

    static native void destroy();

    static native void setRunning(boolean on);

    /** t, dt, nstep, reTau, blown (0/1) */
    static native double[] state();

    /** Nx, Ny, Nz, dx, dz, threads, yc[Ny], yf[Ny+1] */
    static native double[] grid();

    /** plan (Nx*Nz), side (Nx*Ny), end (Nz*Ny) concatenated */
    static native float[] slices(int field, int i, int j, int k);

    /** lo, hi */
    static native double[] range(int field);

    /** averaged (0/1), T, t0, reTau, then yp, Up, Wp, urms, vrms, wrms, uv (Ny each) */
    static native double[] profiles();

    /** n, then t[n], reTau[n], tke[n] */
    static native double[] history();

    static native void startAveragingNow();
}
