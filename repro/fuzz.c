/* Randomized stress fuzzer: mixes normal, duplicate, collinear, huge,
   tiny, zero and (optionally) non-finite coordinates across all cp types,
   both open and closed. Fails on any non-finite bezctx output or sanitizer
   trap. Deterministic via a fixed seed. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spiroentrypoints.h"

static int g_bad;
static int chk(const char *w, double v) {
    if (!isfinite(v)) { g_bad++; if (g_bad <= 10) fprintf(stderr,"  NF %s=%g\n",w,v); return 0; }
    return 1;
}
static void f_mv(bezctx *z,double x,double y,int o){chk("mx",x);chk("my",y);}
static void f_ln(bezctx *z,double x,double y){chk("lx",x);chk("ly",y);}
static void f_qd(bezctx *z,double a,double b,double c,double d)
    {chk("q1",a);chk("q2",b);chk("q3",c);chk("q4",d);}
static void f_cv(bezctx *z,double a,double b,double c,double d,double e,double f)
    {chk("c1",a);chk("c2",b);chk("c3",c);chk("c4",d);chk("c5",e);chk("c6",f);}
static void f_kn(bezctx *z,int k){}

static unsigned long rng;
static double rnd(void) {
    rng = rng * 6364136223846793005UL + 1442695040888963407UL;
    return (double)((rng >> 11) & ((1UL<<53)-1)) / (double)(1UL<<53);
}

static double coord(void) {
    int kind = (int)(rnd() * 10.0);
    double base = rnd() * 200.0 - 100.0;
    switch (kind) {
        case 0: return 0.0;
        case 1: return base * 1e-300;
        case 2: return base * 1e-10;
        case 3: return base;
        case 4: return base * 1e100;
        case 5: return base * 1e200;
        case 6: return base * 1e300;
        case 7: return (rnd() < 0.5 ? 1 : -1) * 1.7976931348623157e308;
        case 8: return rnd() < 0.5 ? INFINITY : -INFINITY;
        default: return NAN;
    }
}

int main(int argc, char **argv) {
    long iters = argc > 1 ? strtol(argv[1], NULL, 10) : 20000;
    rng = 0x12345678abcdefUL;
    static const char closed_ty[] = {'v','o','c'};
    bezctx bc = {f_mv,f_ln,f_qd,f_cv,f_kn};
    long bad_total = 0, ok = 0, rejected = 0;

    for (long it = 0; it < iters; it++) {
        spiro_cp p[10];
        int closed = rnd() < 0.5;
        int n = 3 + (int)(rnd() * 5); /* 3..7 knots */
        if (n > 9) n = 9;
        g_bad = 0;

        if (closed) {
            for (int k = 0; k < n; k++) {
                p[k].x = coord();
                p[k].y = coord();
                p[k].ty = (k == 0) ? 'v' : closed_ty[(int)(rnd()*3)];
            }
            /* sometimes duplicate a knot */
            if (rnd() < 0.3 && n > 2) {
                int d = 1 + (int)(rnd()*(n-1));
                p[d].x = p[d-1].x; p[d].y = p[d-1].y;
            }
            p[n].x = 0; p[n].y = 0; p[n].ty = 'z';
            n++;
        } else {
            p[0].ty = '['; p[n-1].ty = ']';
            for (int k = 0; k < n; k++) {
                p[k].x = coord(); p[k].y = coord();
                if (k > 0 && k < n-1)
                    p[k].ty = closed_ty[(int)(rnd()*3)];
            }
            if (rnd() < 0.3 && n > 2) {
                p[1].x = p[0].x; p[1].y = p[0].y;
            }
        }

        int ncq_choices[] = {0, SPIRO_RETRO_VER1, SPIRO_ARC_MAYBE,
                             SPIRO_QUAD0_TO_BEZIER};
        int ncq = ncq_choices[(int)(rnd()*4)];
        int r = SpiroCPsToBezier2(p, n, ncq, closed, &bc);
        if (g_bad) {
            bad_total += g_bad;
            if (bad_total < 50)
                fprintf(stderr, "iter %ld ret=%d closed=%d n=%d ncq=0x%x\n",
                        it, r, closed, n, ncq);
        }
        if (r == 1) ok++; else rejected++;
    }
    printf("FUZZ iters=%ld ok=%ld clean_reject=%ld nonfinite_outputs=%ld\n",
           iters, ok, rejected, bad_total);
    return bad_total ? 1 : 0;
}
