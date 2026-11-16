/* Benchmark + point-by-point comparison harness for libspiro.
 *
 * Modes:
 *   ./bench dump          print every emitted point (%.17g) for all tiers
 *   ./bench hash          per-contour FNV-1a hashes of all emitted points
 *   ./bench time [reps]   measure solve / emit / total per contour
 *   ./bench batch [reps]  serial vs contour-parallel frame time (pthreads)
 *   ./bench cache [reps]  re-solve-per-frame vs solve-once + emit-per-frame
 *
 * Three outline tiers (control points per closed contour):
 *   S = 24, M = 96, L = 384
 * A "frame" simulates a zoom-preview redraw of one glyph = 8 contours.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>

#include "spiroentrypoints.h"
#include "bezctx.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static unsigned long long rng_state;
static double frand(void)
{
    rng_state ^= rng_state >> 12;
    rng_state ^= rng_state << 25;
    rng_state ^= rng_state >> 27;
    return (double)((rng_state * 0x2545F4914F6CDD1Dull) >> 11) *
	   (1.0 / 9007199254740992.0);
}

/* Deterministic glyph-like closed contour: mostly G4 ('o') points with
   some G2 ('c') and corner ('v') points, on a wobbly rounded shape. */
static spiro_cp *gen_contour(int n, double cx, double cy, double r,
			     unsigned long long seed)
{
    spiro_cp *cp = (spiro_cp *)malloc((size_t)n * sizeof(spiro_cp));
    double ph[4], am[4];
    int i;

    rng_state = seed;
    for (i = 0; i < 4; i++) {
	ph[i] = frand() * 2.0 * M_PI;
	am[i] = 0.10 / (double)(i + 1);
    }
    for (i = 0; i < n; i++) {
	double a = 2.0 * M_PI * (double)i / (double)n;
	double w = 0.0;
	int k;
	for (k = 0; k < 4; k++)
	    w += am[k] * sin((double)(k + 2) * a + ph[k]);
	double rr = r * (1.0 + w);
	cp[i].x = cx + rr * cos(a);
	cp[i].y = cy + rr * 0.82 * sin(a);
	if (i % 48 == 7)
	    cp[i].ty = 'v';
	else if (i % 12 == 3)
	    cp[i].ty = 'c';
	else
	    cp[i].ty = 'o';
    }
    return cp;
}

/* ---- dump context: records every emitted coordinate ---- */
static void d_moveto(bezctx *bc, double x, double y, int is_open)
{
    printf("moveto %.17g %.17g %d\n", x, y, is_open);
}
static void d_lineto(bezctx *bc, double x, double y)
{
    printf("lineto %.17g %.17g\n", x, y);
}
static void d_quadto(bezctx *bc, double x1, double y1, double x2, double y2)
{
    printf("quadto %.17g %.17g %.17g %.17g\n", x1, y1, x2, y2);
}
static void d_curveto(bezctx *bc, double x1, double y1, double x2, double y2,
		      double x3, double y3)
{
    printf("curveto %.17g %.17g %.17g %.17g %.17g %.17g\n",
	   x1, y1, x2, y2, x3, y3);
}
static void d_mark_knot(bezctx *bc, int knot_idx)
{
    printf("mark_knot %d\n", knot_idx);
}

static bezctx dump_ctx = {
    d_moveto, d_lineto, d_quadto, d_curveto, d_mark_knot
};

/* ---- hashing context: FNV-1a over op tags + coordinate bit patterns.
     Doubles as the timing sink and the point-by-point fingerprint. ---- */
typedef struct {
    bezctx base;
    unsigned long long h;
    unsigned long nops;
} hash_ctx;

static void h_mix(hash_ctx *hc, char tag, double v)
{
    unsigned long long b;
    hc->h ^= (unsigned char)tag;
    hc->h *= 1099511628211ull;
    memcpy(&b, &v, 8);
    hc->h ^= b;
    hc->h *= 1099511628211ull;
    hc->nops++;
}
static void h_moveto(bezctx *bc, double x, double y, int is_open)
{
    hash_ctx *hc = (hash_ctx *)bc;
    h_mix(hc, 'm', x); h_mix(hc, (char)is_open, y);
}
static void h_lineto(bezctx *bc, double x, double y)
{
    hash_ctx *hc = (hash_ctx *)bc;
    h_mix(hc, 'l', x); h_mix(hc, 0, y);
}
static void h_quadto(bezctx *bc, double x1, double y1, double x2, double y2)
{
    hash_ctx *hc = (hash_ctx *)bc;
    h_mix(hc, 'q', x1); h_mix(hc, 0, y1);
    h_mix(hc, 0, x2); h_mix(hc, 0, y2);
}
static void h_curveto(bezctx *bc, double x1, double y1, double x2, double y2,
		      double x3, double y3)
{
    hash_ctx *hc = (hash_ctx *)bc;
    h_mix(hc, 'c', x1); h_mix(hc, 0, y1);
    h_mix(hc, 0, x2); h_mix(hc, 0, y2);
    h_mix(hc, 0, x3); h_mix(hc, 0, y3);
}
static void h_mark_knot(bezctx *bc, int knot_idx)
{
    hash_ctx *hc = (hash_ctx *)bc;
    h_mix(hc, 'k', (double)knot_idx);
}

static void hash_ctx_init(hash_ctx *hc)
{
    hc->base.moveto = h_moveto;
    hc->base.lineto = h_lineto;
    hc->base.quadto = h_quadto;
    hc->base.curveto = h_curveto;
    hc->base.mark_knot = h_mark_knot;
    hc->h = 1469598103934665603ull;
    hc->nops = 0;
}

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static const int tier_n[3] = { 24, 96, 384 };
static const char tier_name[3] = { 'S', 'M', 'L' };
#define FRAME_CONTOURS 8
#define NTHREADS 4

static spiro_cp *tier_contour(int t, int c)
{
    return gen_contour(tier_n[t], 500.0 * c, 300.0 * t, 100.0 + 10.0 * c,
		       0x1234ull + 977ull * (unsigned)(t * 64 + c));
}

static void dump_all(void)
{
    int t, c, i;

    for (t = 0; t < 3; t++) {
	int n = tier_n[t];
	for (c = 0; c < FRAME_CONTOURS; c++) {
	    spiro_cp *cp = tier_contour(t, c);
	    spiro_seg *s;
	    printf("tier %c contour %d n %d\n", tier_name[t], c, n);
	    s = run_spiro(cp, n);
	    if (s == NULL) {
		printf("SOLVE FAILED\n");
	    } else {
		for (i = 0; i < n; i++)
		    printf("ks %d %.17g %.17g %.17g %.17g\n", i,
			   s[i].ks[0], s[i].ks[1], s[i].ks[2], s[i].ks[3]);
		free_spiro(s);
	    }
	    if (!SpiroCPsToBezier2(cp, n, 0, 1, &dump_ctx))
		printf("CONVERT FAILED\n");
	    free(cp);
	}
    }
}

static void hash_all(void)
{
    /* serial, batch (contour-parallel) and cached paths must produce
       identical per-contour fingerprints */
    int t, c, w;
    for (t = 0; t < 3; t++) {
	for (c = 0; c < FRAME_CONTOURS; c++) {
	    spiro_cp *cp = tier_contour(t, c);
	    int n = tier_n[t];
	    for (w = 0; w < 3; w++) {
		hash_ctx hc;
		hash_ctx_init(&hc);
		if (w == 0) {
		    SpiroCPsToBezier2(cp, n, 0, 1, &hc.base);
		} else if (w == 1) {
		    /* same call the batch workers make */
		    SpiroCPsToBezier2(cp, n, 0, 1, &hc.base);
		} else {
		    /* cached path: run_spiro0 once, keep s and dm, */
		    /* then spiro_to_bpath0 per frame               */
		    double dm[6];
		    spiro_seg *s = run_spiro0(cp, dm, 0, n);
		    if (s) {
			spiro_to_bpath0(cp, s, dm, 0, n, &hc.base);
			free_spiro(s);
		    }
		}
		printf("tier %c contour %d way %d hash %016llx ops %lu\n",
		       tier_name[t], c, w, hc.h, hc.nops);
	    }
	    free(cp);
	}
    }
}

static double frame_serial(spiro_cp **cps, hash_ctx *hcs, int n)
{
    double t0 = now_sec();
    int c;
    for (c = 0; c < FRAME_CONTOURS; c++)
	SpiroCPsToBezier2(cps[c], n, 0, 1, &hcs[c].base);
    return now_sec() - t0;
}

typedef struct {
    spiro_cp **cps;
    hash_ctx *hcs;
    int n, lo, hi;
} worker_arg;

static void *worker(void *p)
{
    worker_arg *a = (worker_arg *)p;
    int c;
    for (c = a->lo; c < a->hi; c++)
	SpiroCPsToBezier2(a->cps[c], a->n, 0, 1, &a->hcs[c].base);
    return NULL;
}

static double frame_parallel(spiro_cp **cps, hash_ctx *hcs, int n)
{
    pthread_t th[NTHREADS];
    worker_arg args[NTHREADS];
    double t0 = now_sec();
    int t;
    for (t = 0; t < NTHREADS; t++) {
	args[t].cps = cps;
	args[t].hcs = hcs;
	args[t].n = n;
	args[t].lo = t * FRAME_CONTOURS / NTHREADS;
	args[t].hi = (t + 1) * FRAME_CONTOURS / NTHREADS;
	pthread_create(&th[t], NULL, worker, &args[t]);
    }
    for (t = 0; t < NTHREADS; t++)
	pthread_join(th[t], NULL);
    return now_sec() - t0;
}

static double frame_cached(spiro_cp **cps, spiro_seg **segs,
			   double (*dms)[6], hash_ctx *hcs, int n)
{
    double t0 = now_sec();
    int c;
    for (c = 0; c < FRAME_CONTOURS; c++)
	spiro_to_bpath0(cps[c], segs[c], dms[c], 0, n, &hcs[c].base);
    return now_sec() - t0;
}

static void time_all(long reps)
{
    int t, c, r;

    printf("%-4s %6s %12s %12s %12s %12s %12s\n",
	   "tier", "n", "solve_ms", "emit_ms", "frame_ms",
	   "batch4_ms", "cached_ms");
    for (t = 0; t < 3; t++) {
	int n = tier_n[t];
	spiro_cp *cps[FRAME_CONTOURS];
	spiro_seg *segs[FRAME_CONTOURS];
	double dms[FRAME_CONTOURS][6];
	hash_ctx hcs[FRAME_CONTOURS];
	double best_solve = 1e30, best_emit = 1e30, best_full = 1e30;
	double best_batch = 1e30, best_cached = 1e30;

	for (c = 0; c < FRAME_CONTOURS; c++) {
	    cps[c] = tier_contour(t, c);
	    hash_ctx_init(&hcs[c]);
	}
	/* warmup + convergence check */
	for (c = 0; c < FRAME_CONTOURS; c++) {
	    segs[c] = run_spiro0(cps[c], dms[c], 0, n);
	    if (segs[c] == NULL)
		fprintf(stderr, "tier %c contour %d: solve failed\n",
			tier_name[t], c);
	}
	for (r = 0; r < 7; r++) {
	    double dt;
	    long k;

	    /* solve only */
	    dt = 0.;
	    for (k = 0; k <= reps; k++) {
		double t0 = now_sec();
		for (c = 0; c < FRAME_CONTOURS; c++) {
		    double dm[6];
		    spiro_seg *s = run_spiro0(cps[c], dm, 0, n);
		    if (s) free_spiro(s);
		}
		dt += now_sec() - t0;
	    }
	    dt /= (double)(reps + 1);
	    if (dt < best_solve) best_solve = dt;

	    /* emit only (pre-solved segs) */
	    dt = 0.;
	    for (k = 0; k <= reps * 8; k++)
		dt += frame_cached(cps, segs, dms, hcs, n);
	    dt /= (double)(reps * 8 + 1);
	    if (dt < best_emit) best_emit = dt;

	    /* full frame: solve + emit, serial */
	    dt = 0.;
	    for (k = 0; k <= reps; k++)
		dt += frame_serial(cps, hcs, n);
	    dt /= (double)(reps + 1);
	    if (dt < best_full) best_full = dt;

	    /* full frame: solve + emit, contour-parallel */
	    dt = 0.;
	    for (k = 0; k <= reps; k++)
		dt += frame_parallel(cps, hcs, n);
	    dt /= (double)(reps + 1);
	    if (dt < best_batch) best_batch = dt;

	    /* steady-state preview frame: emit from cached solve */
	    dt = 0.;
	    for (k = 0; k <= reps * 8; k++)
		dt += frame_cached(cps, segs, dms, hcs, n);
	    dt /= (double)(reps * 8 + 1);
	    if (dt < best_cached) best_cached = dt;
	}
	printf("%-4c %6d %12.4f %12.4f %12.4f %12.4f %12.4f\n",
	       tier_name[t], n,
	       best_solve * 1e3 / FRAME_CONTOURS,
	       best_emit * 1e3 / FRAME_CONTOURS,
	       best_full * 1e3,
	       best_batch * 1e3,
	       best_cached * 1e3);
	for (c = 0; c < FRAME_CONTOURS; c++) {
	    free_spiro(segs[c]);
	    free(cps[c]);
	}
    }
}

int main(int argc, char **argv)
{
    long reps = 20;

    if (argc > 2) reps = atol(argv[2]);
    if (argc > 1 && strcmp(argv[1], "dump") == 0)
	dump_all();
    else if (argc > 1 && strcmp(argv[1], "hash") == 0)
	hash_all();
    else
	time_all(reps);
    return 0;
}
