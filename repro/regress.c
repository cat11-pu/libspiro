/* Regression: emit a deterministic text dump of every bezctx coordinate
   for a matrix of valid contours / ncq modes. Dumps from the original and
   fixed library are diffed to prove valid contours are unchanged. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spiroentrypoints.h"

static void p_moveto(bezctx *z, double x, double y, int o)
    { printf("M %.17g %.17g %d\n", x, y, o); }
static void p_lineto(bezctx *z, double x, double y)
    { printf("L %.17g %.17g\n", x, y); }
static void p_quadto(bezctx *z, double a, double b, double c, double d)
    { printf("Q %.17g %.17g %.17g %.17g\n", a, b, c, d); }
static void p_curveto(bezctx *z, double a, double b, double c, double d,
                      double e, double f)
    { printf("C %.17g %.17g %.17g %.17g %.17g %.17g\n", a, b, c, d, e, f); }
static void p_knot(bezctx *z, int k) { printf("K %d\n", k); }

static bezctx BC = { p_moveto, p_lineto, p_quadto, p_curveto, p_knot };

static spiro_cp square_closed[] = {
    {0,0,'v'}, {100,0,'o'}, {100,100,'o'}, {0,100,'o'}, {0,0,'z'}
};
static spiro_cp g4_cycle[] = {
    {-100,0,'o'}, {0,100,'o'}, {100,0,'o'}, {0,-100,'o'}
};
static spiro_cp open_o[] = {
    {80,738,'{'}, {749,540,'o'}, {671,309,'o'}, {521,396,'o'},
    {377,333,'o'}, {467,231,'}'}
};
static spiro_cp open_bracket[] = {
    {0,0,'{'}, {100,100,'c'}, {200,200,'['}, {300,200,']'},
    {400,150,'c'}, {300,100,'['}, {200,100,']'}, {150,50,'c'},
    {100,0,'['}, {0,-100,']'}, {-50,-200,'c'}, {-80,-250,'}'}
};
static spiro_cp all_corner[] = {
    {0,0,'v'}, {100,0,'v'}, {100,100,'v'}, {0,100,'v'}, {0,0,'z'}
};
static spiro_cp big_g4[] = {
    {-100e10,0,'o'}, {0,100e10,'o'}, {100e10,0,'o'}, {0,-100e10,'o'}
};
static spiro_cp tiny_g4[] = {
    {-100e-42,0,'o'}, {0,100e-42,'o'}, {100e-42,0,'o'}, {0,-100e-42,'o'}
};

struct { const char *name; spiro_cp *p; int n; int closed; } T[] = {
    {"square_closed", square_closed, 5, 1},
    {"g4_cycle", g4_cycle, 4, 1},
    {"open_o", open_o, 6, 0},
    {"open_bracket", open_bracket, 12, 0},
    {"all_corner", all_corner, 5, 1},
    {"big_g4", big_g4, 4, 1},
    {"tiny_g4", tiny_g4, 4, 1},
};

static const int NCQ[] = {
    0,
    SPIRO_RETRO_VER1,
    SPIRO_CUBIC_MIN_MAYBE,
    SPIRO_ARC_MAYBE,
    SPIRO_ARC_MIN_MAYBE,
    SPIRO_QUAD0_TO_BEZIER,
    SPIRO_INCLUDE_LAST_KNOT,
};

int main(void) {
    for (unsigned t = 0; t < sizeof(T)/sizeof(T[0]); t++) {
        for (unsigned q = 0; q < sizeof(NCQ)/sizeof(NCQ[0]); q++) {
            spiro_cp work[32];
            memcpy(work, T[t].p, T[t].n * sizeof(spiro_cp));
            printf("== %s ncq=0x%x ==\n", T[t].name, NCQ[q]);
            int r = SpiroCPsToBezier2(work, T[t].n, NCQ[q], T[t].closed, &BC);
            printf("RET %d\n", r);
        }
    }
    return 0;
}
