/* Batch conversion, error bound, consistency and boundary tests. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "spiroentrypoints.h"

static spiro_cp path_closed[] = {
    {334, 117, 'v'}, {305, 176, 'v'}, {212, 142, 'c'}, {159, 171, 'c'},
    {224, 237, 'c'}, {347, 335, 'c'}, {202, 467, 'c'}, {81, 429, 'v'},
    {114, 368, 'v'}, {201, 402, 'c'}, {276, 369, 'c'}, {218, 308, 'c'},
    {91, 211, 'c'}, {124, 111, 'c'}, {229, 82, 'c'}, {0, 0, 'z'}
};

static spiro_cp path_open[] = {
    {0, 0, '{'}, {100, 100, 'c'}, {200, 200, 'a'}, {300, 200, 'h'},
    {300, 150, 'c'}, {200, 100, 'a'}, {100, 100, 'h'}, {150, 50, 'c'},
    {100, 0, 'a'}, {0, -100, 'h'}, {50, -100, 'c'}, {20, -150, '}'}
};

static spiro_cp path_circle[] = {
    {-100, 0, 'o'}, {0, 100, 'o'}, {100, 0, 'o'}, {0, -100, 'o'}
};

static spiro_cp path_two[] = {
    {0, 0, 'v'}, {100, 0, 'v'}
};

static int
check_empty(const spiro_batch_result *result)
{
    return result->segments == NULL && result->segment_count == 0 &&
           result->contour_offsets == NULL &&
           result->contour_errors == NULL && result->contour_count == 0;
}

static int
check_geometry(const spiro_batch_result *result, double tolerance)
{
    int c;

    for (c = 0; c < result->contour_count; c++) {
        int begin = result->contour_offsets[c];
        int end = result->contour_offsets[c + 1];
        int i;
        double max_err = 0;
        double cx = 0, cy = 0;
        int have_move = 0;

        if (begin < 0 || end < begin || end > result->segment_count)
            return 0;
        for (i = begin; i < end; i++) {
            const spiro_batch_segment *seg = &result->segments[i];

            if (!isfinite(seg->err) || seg->err < 0 || seg->err > tolerance) {
                printf("bad err c=%d i=%d ty=%c err=%g tol=%g\n",
                       c, i, seg->ty, seg->err, tolerance);
                return 0;
            }
            if (seg->ty == 'm') {
                cx = seg->x1;
                cy = seg->y1;
                have_move = 1;
            } else if (seg->ty == 'l' || seg->ty == 'c') {
                if (!have_move ||
                    hypot(seg->x0 - cx, seg->y0 - cy) >
                    1e-9 * (1 + fabs(cx) + fabs(cy))) {
                    printf("bad continuity c=%d i=%d ty=%c x0=%.17g y0=%.17g cx=%.17g cy=%.17g\n",
                           c, i, seg->ty, seg->x0, seg->y0, cx, cy);
                    return 0;
                }
                if (seg->ty == 'l') {
                    cx = seg->x1;
                    cy = seg->y1;
                } else {
                    cx = seg->x3;
                    cy = seg->y3;
                }
                if (seg->err > max_err)
                    max_err = seg->err;
            } else if (seg->ty != 'k') {
                printf("bad type c=%d i=%d ty=%d\n", c, i, seg->ty);
                return 0;
            }
        }
        if (fabs(max_err - result->contour_errors[c]) > 0) {
            printf("bad contour max c=%d got=%.17g expect=%.17g\n",
                   c, max_err, result->contour_errors[c]);
            return 0;
        }
    }
    return 1;
}

static int
check_sequential(const spiro_contour *contours, int ncontours, int ncq,
                 double tolerance)
{
    spiro_batch_result batch;
    int c;

    if (!SpiroCPsToBezierBatch(contours, ncontours, ncq, tolerance, &batch)) {
        printf("FAIL: batch conversion at tolerance %g\n", tolerance);
        return 0;
    }
    if (!check_geometry(&batch, tolerance)) {
        printf("FAIL: batch geometry/error check at tolerance %g\n", tolerance);
        SpiroFreeBatchResult(&batch);
        return 0;
    }

    for (c = 0; c < ncontours; c++) {
        spiro_batch_result one;
        int begin = batch.contour_offsets[c];
        int count = batch.contour_offsets[c + 1] - begin;

        if (!SpiroCPsToBezierTol(contours[c].points, contours[c].n, ncq,
                                 contours[c].is_closed, tolerance, &one)) {
            printf("FAIL: single conversion, contour %d\n", c);
            SpiroFreeBatchResult(&batch);
            return 0;
        }
        if (one.segment_count != count ||
            memcmp(one.segments, batch.segments + begin,
                   (size_t)count * sizeof(*one.segments)) != 0 ||
            one.contour_errors[0] != batch.contour_errors[c]) {
            printf("FAIL: batch != sequential, contour %d, tolerance %g\n",
                   c, tolerance);
            SpiroFreeBatchResult(&one);
            SpiroFreeBatchResult(&batch);
            return 0;
        }
        SpiroFreeBatchResult(&one);
    }

    {
        double all_max = 0;
        int geometric = 0;
        for (c = 0; c < ncontours; c++) {
            if (batch.contour_errors[c] > all_max)
                all_max = batch.contour_errors[c];
        }
        for (c = 0; c < batch.segment_count; c++) {
            if (batch.segments[c].ty == 'l' || batch.segments[c].ty == 'c')
                ++geometric;
        }
        printf("consistency: ncq=0x%x tol=%-8g commands=%-4d "
               "geometry=%-4d max_err=%.9g\n",
               ncq, tolerance, batch.segment_count, geometric, all_max);
    }
    SpiroFreeBatchResult(&batch);
    return 1;
}

static int
test_consistency(void)
{
    static const double tolerances[] = {1.0, 0.1, 0.01, 0.001, 0.0001};
    spiro_contour contours[] = {
        {path_closed, 16, 1},
        {path_open, 12, 0},
        {path_circle, 4, 1}
    };
    int i;

    for (i = 0; i < (int)(sizeof(tolerances) / sizeof(tolerances[0])); i++) {
        if (!check_sequential(contours, 3, 0, tolerances[i]))
            return 0;
        if (!check_sequential(contours, 3, SPIRO_RETRO_VER1, tolerances[i]))
            return 0;
    }
    return 1;
}

static int
test_linear_output(void)
{
    enum { COPIES = 16 };
    spiro_contour many[COPIES];
    spiro_batch_result one, batch;
    int i;

    for (i = 0; i < COPIES; i++) {
        many[i].points = path_closed;
        many[i].n = 16;
        many[i].is_closed = 1;
    }
    if (!SpiroCPsToBezierTol(path_closed, 16, 0, 1, 0.1, &one) ||
        !SpiroCPsToBezierBatch(many, COPIES, 0, 0.1, &batch)) {
        printf("FAIL: linear-output setup\n");
        SpiroFreeBatchResult(&one);
        SpiroFreeBatchResult(&batch);
        return 0;
    }
    if (batch.segment_count != one.segment_count * COPIES) {
        printf("FAIL: repeated contours are not linear\n");
        SpiroFreeBatchResult(&one);
        SpiroFreeBatchResult(&batch);
        return 0;
    }
    for (i = 0; i < COPIES; i++) {
        int begin = batch.contour_offsets[i];
        int count = batch.contour_offsets[i + 1] - begin;
        if (count != one.segment_count ||
            batch.contour_errors[i] != one.contour_errors[0] ||
            memcmp(batch.segments + begin, one.segments,
                   (size_t)count * sizeof(*one.segments)) != 0) {
            printf("FAIL: repeated contour %d differs\n", i);
            SpiroFreeBatchResult(&one);
            SpiroFreeBatchResult(&batch);
            return 0;
        }
    }
    printf("linearity: contours=%-3d points=%-4d segments=%-4d\n",
           COPIES, COPIES * 16, batch.segment_count);
    SpiroFreeBatchResult(&one);
    SpiroFreeBatchResult(&batch);
    return 1;
}

static int
expect_failure(const spiro_contour *contours, int ncontours, int ncq,
               double tolerance, const char *what)
{
    spiro_batch_result result;

    if (SpiroCPsToBezierBatch(contours, ncontours, ncq, tolerance, &result) ||
        !check_empty(&result)) {
        printf("FAIL: expected rejection: %s\n", what);
        SpiroFreeBatchResult(&result);
        return 0;
    }
    return 1;
}

static int
test_boundaries(void)
{
    spiro_contour good = {path_closed, 16, 1};
    spiro_contour bad_mid[] = {
        {path_closed, 16, 1},
        {path_open, 1, 0},
        {path_circle, 4, 1}
    };
    spiro_cp bad_point[2] = {{0, 0, 'v'}, {NAN, 1, 'v'}};
    spiro_contour bad_finite = {bad_point, 2, 0};
    spiro_batch_result empty;
    int ncq_bad[] = {
        SPIRO_INTERNAL_BEZCTX,
        SPIRO_REVERSE_SRC,
        SPIRO_QUAD0_TO_BEZIER
    };
    int i;

    if (!SpiroCPsToBezierBatch(NULL, 0, 0, 1.0, &empty) ||
        !check_empty(&empty) || empty.segment_count != 0) {
        printf("FAIL: zero-contour batch\n");
        SpiroFreeBatchResult(&empty);
        return 0;
    }
    SpiroFreeBatchResult(&empty);
    if (SpiroCPsToBezierBatch(NULL, 0, 0, 1.0, NULL) ||
        SpiroCPsToBezierBatch(NULL, 1, 0, 1.0, &empty) ||
        !check_empty(&empty)) {
        printf("FAIL: NULL argument handling\n");
        SpiroFreeBatchResult(&empty);
        return 0;
    }

    if (!expect_failure(&good, 1, 0, 0.0, "zero tolerance") ||
        !expect_failure(&good, 1, 0, -1.0, "negative tolerance") ||
        !expect_failure(&good, 1, 0, NAN, "NaN tolerance") ||
        !expect_failure(&good, 1, 0, INFINITY, "infinite tolerance"))
        return 0;
    for (i = 0; i < (int)(sizeof(ncq_bad) / sizeof(ncq_bad[0])); i++) {
        if (!expect_failure(&good, 1, ncq_bad[i], 1.0, "unsupported ncq"))
            return 0;
    }
    if (!expect_failure(bad_mid, 3, 0, 1.0, "invalid middle contour") ||
        !expect_failure(&bad_finite, 1, 0, 1.0, "non-finite point"))
        return 0;

    {
        spiro_contour one_point = {path_two, 1, 0};
        if (!expect_failure(&one_point, 1, 0, 1.0, "one-point contour"))
            return 0;
    }
    printf("boundaries: empty batch, NULL, tolerance, ncq, atomic failure\n");
    return 1;
}

static int
test_last_knot(void)
{
    spiro_batch_result result;
    int i, knots = 0;

    if (!SpiroCPsToBezierTol(path_open, 12,
                             SPIRO_INCLUDE_LAST_KNOT, 0, 0.1, &result)) {
        printf("FAIL: last-knot conversion\n");
        return 0;
    }
    for (i = 0; i < result.segment_count; i++) {
        if (result.segments[i].ty == 'k')
            ++knots;
    }
    if (knots != 9) {
        printf("FAIL: last-knot count %d\n", knots);
        SpiroFreeBatchResult(&result);
        return 0;
    }
    SpiroFreeBatchResult(&result);
    return 1;
}

int
main(void)
{
    if (!test_boundaries() || !test_last_knot() ||
        !test_consistency() || !test_linear_output())
        return 1;
    printf("pass: batch conversion\n");
    return 0;
}
