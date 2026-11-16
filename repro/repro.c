/* Reproduction harness for degenerate / duplicate / huge-coordinate inputs.
 * Tracks every bezctx output coordinate and flags NaN/Inf, plus crashes.
 * Build: see build.sh
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spiroentrypoints.h"

static int g_nan_count;
static int g_call_count;

static int check(const char *what, double v) {
    if (!isfinite(v)) {
        g_nan_count++;
        fprintf(stderr, "  [NON-FINITE] %s = %g\n", what, v);
        return 0;
    }
    return 1;
}

static void cap_moveto(bezctx *bc, double x, double y, int is_open) {
    g_call_count++;
    check("moveto.x", x); check("moveto.y", y);
}
static void cap_lineto(bezctx *bc, double x, double y) {
    g_call_count++;
    check("lineto.x", x); check("lineto.y", y);
}
static void cap_quadto(bezctx *bc, double x1, double y1, double x2, double y2) {
    g_call_count++;
    check("quadto.x1", x1); check("quadto.y1", y1);
    check("quadto.x2", x2); check("quadto.y2", y2);
}
static void cap_curveto(bezctx *bc, double x1, double y1, double x2, double y2,
                        double x3, double y3) {
    g_call_count++;
    check("curveto.x1", x1); check("curveto.y1", y1);
    check("curveto.x2", x2); check("curveto.y2", y2);
    check("curveto.x3", x3); check("curveto.y3", y3);
}
static void cap_mark_knot(bezctx *bc, int knot_idx) {}

static bezctx *new_cap(void) {
    bezctx *bc = calloc(1, sizeof(bezctx));
    bc->moveto = cap_moveto;
    bc->lineto = cap_lineto;
    bc->quadto = cap_quadto;
    bc->curveto = cap_curveto;
    bc->mark_knot = cap_mark_knot;
    return bc;
}

typedef struct {
    const char *name;
    int closed;
    spiro_cp pts[12];
    int n;
} test_case;

#define CP(x,y,t) {(x),(y),(t)}
#define ENDZ CP(0,0,'z')

static test_case cases[] = {
    /* baseline: a normal valid closed contour */
    {"valid_closed_square", 1,
     { CP(0,0,'v'), CP(100,0,'o'), CP(100,100,'o'),
       CP(0,100,'o'), ENDZ }, 5},
    /* baseline: normal open contour */
    {"valid_open", 0,
     { CP(0,0,'['), CP(50,100,'c'), CP(100,0,']') }, 3},

    /* duplicate consecutive points (zero-length chords) */
    {"dup_closed_all_same", 1,
     { CP(10,10,'v'), CP(10,10,'o'), CP(10,10,'o'),
       CP(10,10,'o'), ENDZ }, 5},
    {"dup_closed_two_same", 1,
     { CP(0,0,'v'), CP(100,0,'o'), CP(100,0,'o'),
       CP(0,100,'o'), ENDZ }, 5},
    {"dup_open_same", 0,
     { CP(5,5,'['), CP(5,5,'c'), CP(5,5,']') }, 3},

    /* fully degenerate: 2-point closed */
    {"two_point_closed", 1,
     { CP(0,0,'v'), CP(10,10,'v'), ENDZ }, 3},

    /* huge coordinates */
    {"huge_closed_1e150", 1,
     { CP(-1e150,-1e150,'v'), CP(1e150,-1e150,'o'),
       CP(1e150,1e150,'o'), CP(-1e150,1e150,'o'), ENDZ }, 5},
    {"huge_closed_1e308", 1,
     { CP(-1e308,-1e308,'v'), CP(1e308,-1e308,'o'),
       CP(1e308,1e308,'o'), CP(-1e308,1e308,'o'), ENDZ }, 5},
    {"huge_dup_1e200", 1,
     { CP(1e200,1e200,'v'), CP(1e200,1e200,'o'),
       CP(1e200,1e200,'o'), CP(1e200,1e200,'o'), ENDZ }, 5},

    /* explicit non-finite inputs */
    {"inf_coords", 1,
     { CP(0,0,'v'), CP(INFINITY,0,'o'),
       CP(100,100,'o'), CP(0,100,'o'), ENDZ }, 5},
    {"nan_coords", 1,
     { CP(0,0,'v'), CP(NAN,0,'o'),
       CP(100,100,'o'), CP(0,100,'o'), ENDZ }, 5},
};

int main(void) {
    int n_cases = (int)(sizeof(cases) / sizeof(cases[0]));
    int total_nan = 0, total_calls = 0, total_fail_ret = 0;

    for (int c = 0; c < n_cases; c++) {
        test_case *tc = &cases[c];
        bezctx *bc = new_cap();
        g_nan_count = 0;
        g_call_count = 0;

        printf("CASE %-22s n=%d closed=%d\n", tc->name, tc->n, tc->closed);
        int ret = SpiroCPsToBezier2(tc->pts, tc->n, 0, tc->closed, bc);
        printf("  return=%d  bezctx_calls=%d  non_finite_outputs=%d\n",
               ret, g_call_count, g_nan_count);

        total_nan += g_nan_count;
        total_calls += g_call_count;
        if (ret != 1 && g_nan_count == 0)
            total_fail_ret++;
        free(bc);
    }
    printf("SUMMARY non_finite=%d total_bezctx_calls=%d clean_fail_cases=%d\n",
           total_nan, total_calls, total_fail_ret);
    return total_nan ? 2 : 0;
}
