/*
libspiro - A sharable library of Spiro formula and functions.

This file is a fork from ppedit for use by the libspiro project.
Please see Changelog or git history for description of changes.
=============================================================
Copyright (C) 2007... Raph Levien (ppedit)
Copyright (C) 2007... George Williams (libspiro fork)
Copyright (C) 2013... Joe Da Silva (improvements)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
*/
#include "spiroentrypoints.h"
#include "spiro_intf.h"
#include "spiro-config.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* These six functions are kept for backwards compatibility */
LS_DLL_EXPORT
void SpiroCPsToBezier(spiro_cp *spiros,int n,int isclosed,bezctx *bc) {
    SpiroCPsToBezier2(spiros,n,SPIRO_RETRO_VER1,isclosed,bc);
}
LS_DLL_EXPORT
void TaggedSpiroCPsToBezier(spiro_cp *spiros,bezctx *bc) {
    TaggedSpiroCPsToBezier2(spiros,SPIRO_RETRO_VER1,bc);
}
LS_DLL_EXPORT
int SpiroCPsToBezier0(spiro_cp *spiros,int n,int isclosed,bezctx *bc) {
    return SpiroCPsToBezier2(spiros,n,SPIRO_RETRO_VER1,isclosed,bc);
}
LS_DLL_EXPORT
int TaggedSpiroCPsToBezier0(spiro_cp *spiros,bezctx *bc) {
    return TaggedSpiroCPsToBezier2(spiros,SPIRO_RETRO_VER1,bc);
}
LS_DLL_EXPORT
void SpiroCPsToBezier1(spiro_cp *spiros,int n,int isclosed,bezctx *bc,int *done) {
    *done = SpiroCPsToBezier2(spiros,n,SPIRO_RETRO_VER1,isclosed,bc);
}
LS_DLL_EXPORT
void TaggedSpiroCPsToBezier1(spiro_cp *spiros,bezctx *bc,int *done) {
    *done = TaggedSpiroCPsToBezier2(spiros,SPIRO_RETRO_VER1,bc);
}


LS_DLL_EXPORT
int SpiroCPsToBezier2(spiro_cp *spiros,int n,int ncq,int isclosed,bezctx *bc) {
    double dm[6];
    spiro_seg *s;

    if ( n<=0 || ncq<0 ) return 0; /* invalid inputs */
    if ( n > 1 && (ncq & SPIRO_REVERSE_SRC) && spiroreverse(spiros,n) )
	return 0;

    if ( isclosed )
	s = run_spiro0(spiros,dm,ncq,n);
    else {
	char oldty_start = spiros[0].ty;
	char oldty_end   = spiros[n-1].ty;
	spiros[0].ty = '{';
	spiros[n-1].ty = '}';
	s = run_spiro0(spiros,dm,ncq,n);
	spiros[n-1].ty = oldty_end;
	spiros[0].ty = oldty_start;
    }
    if ( s ) {
	if ( spiro_to_bpath0(spiros,s,dm,ncq,n,bc) ) {
	    free_spiro(s);
	    return 1; /* success */
	}
	free_spiro(s);
    }
    return 0; /* spiro did not converge or encountered non-finite values */
}

LS_DLL_EXPORT
int TaggedSpiroCPsToBezier2(spiro_cp *spiros,int ncq,bezctx *bc) {
    double dm[6];
    spiro_seg *s;
    int n;

    for ( n=0; spiros[n].ty!='z' && spiros[n].ty!='}'; ++n );
    if ( spiros[n].ty == '}' ) ++n;

    if ( n<=0 || ncq<0 ) return 0; /* invalid inputs */
    if ( n > 1 && (ncq & SPIRO_REVERSE_SRC) && spiroreverse(spiros,n) )
	return 0;
    s = run_spiro0(spiros,dm,ncq,n);
    if ( s ) {
	if ( spiro_to_bpath0(spiros,s,dm,ncq,n,bc) ) {
	    free_spiro(s);
	    return 1; /* success */
	}
	free_spiro(s);
    }
    return 0; /* spiro did not converge or encountered non-finite values */
}

static int
spiro_batch_limited_ncq(int ncq)
{
    int supported = SPIRO_INCLUDE_LAST_KNOT | SPIRO_RETRO_VER1;

    return ncq >= 0 && (ncq & ~supported) == 0;
}

static int
spiro_batch_valid_contour(const spiro_contour *contour)
{
    int i;

    if (contour == NULL || contour->points == NULL || contour->n < 2)
        return 0;
    for (i = 0; i < contour->n; i++) {
        if (!IS_FINITE(contour->points[i].x) ||
            !IS_FINITE(contour->points[i].y))
            return 0;
    }
    return 1;
}

typedef struct {
    spiro_batch_result *result;
    double max_err;
} spiro_batch_builder;

static int
spiro_batch_reserve(spiro_batch_result *result, int needed)
{
    int capacity;
    spiro_batch_segment *segments;

    if (needed <= result->segment_capacity)
        return 1;
    capacity = result->segment_capacity > 0 ? result->segment_capacity : 64;
    while (capacity < needed) {
        if (capacity > INT_MAX / 2)
            return 0;
        capacity *= 2;
    }
    segments = (spiro_batch_segment *)realloc(
        result->segments, (size_t)capacity * sizeof(*segments));
    if (segments == NULL)
        return 0;
    result->segments = segments;
    result->segment_capacity = capacity;
    return 1;
}

static int
spiro_batch_emit(void *ctx, const spiro_batch_segment *segment)
{
    spiro_batch_builder *builder = (spiro_batch_builder *)ctx;
    spiro_batch_result *result = builder->result;

    if (result->segment_count == INT_MAX ||
        !IS_FINITE(segment->err) || segment->err < 0)
        return 0;
    if (!spiro_batch_reserve(result, result->segment_count + 1))
        return 0;
    result->segments[result->segment_count++] = *segment;
    if (segment->ty == 'l' || segment->ty == 'c') {
        if (segment->err > builder->max_err)
            builder->max_err = segment->err;
    }
    return 1;
}

static int
spiro_batch_convert_one(const spiro_contour *contour, int ncq,
                        double tolerance, spiro_batch_result *result,
                        double *max_err)
{
    spiro_cp *copy = NULL;
    const spiro_cp *src = contour->points;
    spiro_seg *segs;
    spiro_batch_builder builder;
    double dm[6];
    int ok;

    if (!contour->is_closed) {
        copy = (spiro_cp *)malloc((size_t)contour->n * sizeof(*copy));
        if (copy == NULL)
            return 0;
        memcpy(copy, contour->points, (size_t)contour->n * sizeof(*copy));
        copy[0].ty = SPIRO_OPEN_CONTOUR;
        copy[contour->n - 1].ty = SPIRO_END_OPEN_CONTOUR;
        src = copy;
    }

    segs = run_spiro0(src, dm, ncq, contour->n);
    if (segs == NULL) {
        free(copy);
        return 0;
    }

    builder.result = result;
    builder.max_err = 0;
    ok = spiro_to_bpath_tol0(contour->points, segs, dm, ncq, contour->n,
                             tolerance, spiro_batch_emit, &builder);
    free_spiro(segs);
    free(copy);
    if (!ok)
        return 0;
    *max_err = builder.max_err;
    return 1;
}

LS_DLL_EXPORT
int
SpiroCPsToBezierTol(const spiro_cp *points, int n, int ncq,
                    int is_closed, double tolerance,
                    spiro_batch_result *result)
{
    spiro_contour contour;

    contour.points = points;
    contour.n = n;
    contour.is_closed = is_closed;
    return SpiroCPsToBezierBatch(&contour, 1, ncq, tolerance, result);
}

LS_DLL_EXPORT
int
SpiroCPsToBezierBatch(const spiro_contour *contours, int ncontours,
                      int ncq, double tolerance,
                      spiro_batch_result *result)
{
    int i;

    if (result == NULL)
        return 0;
    memset(result, 0, sizeof(*result));
    if (ncontours < 0 || !(tolerance > 0) || !IS_FINITE(tolerance) ||
        !spiro_batch_limited_ncq(ncq))
        return 0;
    if (ncontours == 0)
        return 1;
    if (contours == NULL)
        return 0;
    for (i = 0; i < ncontours; i++) {
        if (!spiro_batch_valid_contour(&contours[i]))
            return 0;
    }

    result->contour_offsets = (int *)calloc((size_t)ncontours + 1,
                                            sizeof(*result->contour_offsets));
    result->contour_errors = (double *)calloc((size_t)ncontours,
                                              sizeof(*result->contour_errors));
    if (result->contour_offsets == NULL || result->contour_errors == NULL) {
        SpiroFreeBatchResult(result);
        return 0;
    }
    result->contour_count = ncontours;

    for (i = 0; i < ncontours; i++) {
        result->contour_offsets[i] = result->segment_count;
        if (!spiro_batch_convert_one(&contours[i], ncq, tolerance, result,
                                     &result->contour_errors[i])) {
            SpiroFreeBatchResult(result);
            return 0;
        }
    }
    result->contour_offsets[ncontours] = result->segment_count;
    return 1;
}

LS_DLL_EXPORT
void
SpiroFreeBatchResult(spiro_batch_result *result)
{
    if (result == NULL)
        return;
    free(result->segments);
    free(result->contour_offsets);
    free(result->contour_errors);
    memset(result, 0, sizeof(*result));
}
