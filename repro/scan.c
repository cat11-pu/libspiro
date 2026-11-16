/* Boundary scan: enumerate degenerate / duplicate / huge-coordinate inputs.
 * Reports any case where the library returns success but emits a
 * non-finite bezctx coordinate, or where sanitizers trap.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spiroentrypoints.h"

static int g_bad;
static int g_calls;

static int chk(const char *w, double v) {
    if (!isfinite(v)) {
        g_bad++;
        if (g_bad <= 25)
            fprintf(stderr, "    NON-FINITE %s=%g\n", w, v);
        return 0;
    }
    return 1;
}

static void m_moveto(bezctx *bc, double x, double y, int o) {
    g_calls++; chk("moveto.x", x); chk("moveto.y", y);
}
static void m_lineto(bezctx *bc, double x, double y) {
    g_calls++; chk("lineto.x", x); chk("lineto.y", y);
}
static void m_quadto(bezctx *bc, double a, double b, double c, double d) {
    g_calls++; chk("q.1", a); chk("q.2", b); chk("q.3", c); chk("q.4", d);
}
static void m_curveto(bezctx *bc, double a, double b, double c, double d,
                      double e, double f) {
    g_calls++;
    chk("c.1", a); chk("c.2", b); chk("c.3", c);
    chk("c.4", d); chk("c.5", e); chk("c.6", f);
}
static void m_knot(bezctx *bc, int k) {}

static bezctx *new_bc(void) {
    bezctx *bc = calloc(1, sizeof(bezctx));
    bc->moveto = m_moveto; bc->lineto = m_lineto;
    bc->quadto = m_quadto; bc->curveto = m_curveto;
    bc->mark_knot = m_knot;
    return bc;
}

static const char *TY_CLOSED[] = {"v", "o", "c"};
static const char TY_OPEN[3][3] = { {'[','c',']'}, {'[','o',']'}, {'[','v',']'} };

static const double MAG[] = {
    0.0, 1e-320, 1e-300, 1e-100, 1e-10, 1.0, 1e10, 1e100,
    1e150, 1e200, 1e300, 1e308,
};
#define NMAG ((int)(sizeof(MAG)/sizeof(MAG[0])))

static double val(double mag, int quad, int axis) {
    static const double sx[4] = { 1, -1, 1, -1 };
    static const double sy[4] = { 1, 1, -1, -1 };
    return mag * (axis == 0 ? sx[quad] : sy[quad]);
}

static long total, bad_cases;

static long cat_total[5], cat_bad[5], cat_ok[5], cat_reject[5];

/* category-aware variant that records per-bucket statistics */
static void run_cat(int cat, const char *label, spiro_cp *p, int n, int closed) {
    bezctx *bc = new_bc();
    g_bad = 0; g_calls = 0;
    cat_total[cat]++; total++;
    int r = SpiroCPsToBezier2(p, n, 0, closed, bc);
    if (g_bad) { cat_bad[cat]++; bad_cases++;
        fprintf(stderr, "BAD label=%s ret=%d bad=%d\n", label, r, g_bad); }
    if (r == 1) cat_ok[cat]++; else cat_reject[cat]++;
    free(bc);
}

int main(void) {
    spiro_cp p[12];
    const char *catname[5] = {
        "A:closed-magnitude", "B:duplicate-chord", "C:open-magnitude",
        "D:nonfinite-inject", "TOTAL"};

    /* ---- scan A: closed 4-cp contours, coordinate magnitude sweep ---- */
    for (int mi = 0; mi < NMAG; mi++) {
        for (int quad = 0; quad < 4; quad++) {
            for (int ti = 0; ti < 3; ti++) {
                double m = MAG[mi];
                char t = TY_CLOSED[ti][0];
                p[0] = (spiro_cp){ val(m,quad,0), val(m,quad,1), 'v' };
                p[1] = (spiro_cp){ val(m,quad,0) + m, val(m,quad,1), t };
                p[2] = (spiro_cp){ val(m,quad,0) + m, val(m,quad,1) + m, t };
                p[3] = (spiro_cp){ val(m,quad,0), val(m,quad,1) + m, t };
                p[4] = (spiro_cp){ 0, 0, 'z' };
                char lab[64];
                snprintf(lab, sizeof lab, "closed mag=%.0e ty=%c q=%d", m, t, quad);
                run_cat(0, lab, p, 5, 1);
            }
        }
    }

    /* ---- scan B: closed contours with one duplicated (zero) chord ---- */
    for (int mi = 0; mi < NMAG; mi++) {
        double m = MAG[mi];
        for (int dup = 0; dup < 4; dup++) {
            for (int ti = 0; ti < 3; ti++) {
                char t = TY_CLOSED[ti][0];
                for (int k = 0; k < 4; k++)
                    p[k] = (spiro_cp){ val(m,0,0), val(m,0,1), t };
                p[0].ty = 'v';
                p[1].x += m;
                p[2].x += m; p[2].y += m;
                p[3].y += m;
                /* force one chord to zero length */
                p[(dup + 1) % 4].x = p[dup].x;
                p[(dup + 1) % 4].y = p[dup].y;
                p[4] = (spiro_cp){ 0, 0, 'z' };
                char lab[64];
                snprintf(lab, sizeof lab, "dup mag=%.0e ty=%c at=%d", m, t, dup);
                run_cat(1, lab, p, 5, 1);
            }
        }
    }

    /* ---- scan C: open 3-cp contours, magnitude + type sweep ---- */
    for (int mi = 0; mi < NMAG; mi++) {
        double m = MAG[mi];
        for (int ti = 0; ti < 3; ti++) {
            for (int mid = 0; mid < 4; mid++) {
                p[0] = (spiro_cp){ val(m,0,0), val(m,0,1), TY_OPEN[ti][0] };
                p[1] = (spiro_cp){ val(m,mid,0), val(m,mid,1), TY_OPEN[ti][1] };
                p[2] = (spiro_cp){ val(m,0,0) + m, val(m,0,1), TY_OPEN[ti][2] };
                char lab[64];
                snprintf(lab, sizeof lab, "open mag=%.0e ty=%c mid=%d", m,
                         TY_OPEN[ti][1], mid);
                run_cat(2, lab, p, 3, 0);
            }
        }
    }

    /* ---- scan D: non-finite coordinate injection ---- */
    double specials[] = { INFINITY, -INFINITY, NAN };
    const char *sname[] = {"+inf", "-inf", "nan"};
    for (int si = 0; si < 3; si++) {
        for (int pos = 0; pos < 4; pos++) {
            for (int axis = 0; axis < 2; axis++) {
                for (int k = 0; k < 4; k++)
                    p[k] = (spiro_cp){ 10.0 * k, 10.0 * (k % 2), 'o' };
                p[0].ty = 'v';
                if (axis == 0) p[pos].x = specials[si];
                else p[pos].y = specials[si];
                p[4] = (spiro_cp){ 0, 0, 'z' };
                char lab[64];
                snprintf(lab, sizeof lab, "special %s pos=%d axis=%d",
                         sname[si], pos, axis);
                run_cat(3, lab, p, 5, 1);
            }
        }
    }

    printf("SCAN total=%ld bad_nonfinite=%ld\n", total, bad_cases);
    printf("\n%-20s %8s %8s %8s %8s\n", "category", "total", "ok", "reject", "nonfinite");
    for (int c = 0; c < 4; c++)
        printf("%-20s %8ld %8ld %8ld %8ld\n", catname[c],
               cat_total[c], cat_ok[c], cat_reject[c], cat_bad[c]);
    long st=0,sok=0,srj=0,sbd=0;
    for (int c=0;c<4;c++){st+=cat_total[c];sok+=cat_ok[c];srj+=cat_reject[c];sbd+=cat_bad[c];}
    printf("%-20s %8ld %8ld %8ld %8ld\n", catname[4], st, sok, srj, sbd);
    return bad_cases ? 1 : 0;
}
