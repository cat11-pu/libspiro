# libspiro curve-solver optimization bench

Harness for the zoom-preview workload: one "frame" = solve + flatten of
one glyph made of 8 closed contours, at three outline tiers
(S = 24, M = 96, L = 384 control points per contour).

## Build

    # baseline (pristine spiro.c from git HEAD)
    git show HEAD:spiro.c > bench/spiro_orig.c
    gcc -O2 -I. -Ibench bench/spiro_orig.c spiroentrypoints.c bezctx.c \
        bench/bench.c -lm -lpthread -o bench_base
    # optimized
    gcc -O2 -I. -Ibench spiro.c spiroentrypoints.c bezctx.c \
        bench/bench.c -lm -lpthread -o bench_opt

(`bench/spiro-config.h` is a minimal stand-in for the autotools-generated
config header; the real build uses `./configure && make`.)

## Run

    ./bench_base dump  > base.txt   # every emitted point, %.17g
    ./bench_opt  dump  > opt.txt
    diff base.txt opt.txt           # must be empty: bit-identical points
    ./bench_opt  hash               # per-contour FNV-1a fingerprints for
                                    # serial / batch / cached paths
    ./bench_opt  time 40            # timings (best of 7 x 41 reps)

## What the columns mean

- `solve_ms`   per-contour `run_spiro0()` (the Newton solve)
- `emit_ms`    per-contour `spiro_to_bpath0()` from a solved `spiro_seg`
- `frame_ms`   8 contours, solve + emit, serial
- `batch4_ms`  same frame, contours split over 4 pthreads (existing
               public API only; libspiro calls are thread-safe, see
               tests/call-testm.c)
- `cached_ms`  steady-state preview frame: solve once, keep the
               `spiro_seg` + `dm`, re-emit per frame

## Caching pattern (zoom preview)

The solve does not depend on the zoom factor, so per-frame re-solving is
unnecessary. Solve once per contour, reuse across frames:

    double dm[6];
    spiro_seg *s = run_spiro0(cps, dm, ncq, n);     /* once */
    /* per frame: */
    spiro_to_bpath0(cps, s, dm, ncq, n, bc);
    /* on control-point edit: free_spiro(s) and re-solve */

Use `run_spiro0`/`spiro_to_bpath0` (not the deprecated `run_spiro`/
`spiro_to_bpath`): the deprecated pair forces the identity transform
(`SPIRO_RETRO_VER1`) and is *not* bit-identical to `SpiroCPsToBezier2`
with `ncq == 0`, which auto-normalizes coordinates via `dm`.
