/* spirobatch.h
libspiro - batch Spiro contour conversion with per-segment error bounds.
*/
#ifndef _SPIROBATCH_H
#define _SPIROBATCH_H

#include "spiro.h"

#ifdef __cplusplus
extern "C" {
#endif

/* One contour in a batch.  Open contours do not need '{' / '}' tags; the
   caller supplies is_closed instead.  Closed contours do not need 'z'. */
typedef struct {
    const spiro_cp *points;
    int n;
    int is_closed;
} spiro_contour;

/* A flattened Bezier command.  Field meanings follow curve_data:
     'm': move to (x1,y1)
     'k': knot marker, no geometry
     'l': line from (x0,y0) to (x1,y1)
     'c': cubic from (x0,y0), controls (x1,y1)/(x2,y2), to (x3,y3)
   err is a conservative upper bound for the distance between this command
   and its source Spiro segment, in the same units as the input coordinates.
   Move and knot commands have err == 0. */
typedef struct {
    char ty;
    double x0, y0;
    double x1, y1;
    double x2, y2;
    double x3, y3;
    double err;
} spiro_batch_segment;

/* Segments are flattened in contour order.  Contour c owns
   segments[contour_offsets[c] .. contour_offsets[c + 1]).
   contour_errors[c] is the maximum geometric err in that range. */
typedef struct {
    spiro_batch_segment *segments;
    int segment_count;
    int segment_capacity;       /* internal; do not modify */
    int *contour_offsets;       /* contour_count + 1 entries */
    double *contour_errors;     /* contour_count entries */
    int contour_count;
} spiro_batch_result;

/* Convert one contour.  Returns 1 on success, 0 on invalid input, failure to
   converge, insufficient memory, or inability to meet tolerance.

   Only cubic output is supported.  ncq may contain RETRO and
   INCLUDE_LAST_KNOT; INTERNAL_BEZCTX, REVERSE_SRC and the arc/quad output
   selectors are rejected.  tolerance must be finite and greater than zero. */
int SpiroCPsToBezierTol(const spiro_cp *points, int n, int ncq,
                        int is_closed, double tolerance,
                        spiro_batch_result *result);

/* Batch form of SpiroCPsToBezierTol().  On success, each contour is converted
   by exactly the same code path as the single-contour entry point.  A failure
   is atomic: result is left empty.  Call SpiroFreeBatchResult() before reusing
   a populated result for another conversion. */
int SpiroCPsToBezierBatch(const spiro_contour *contours, int ncontours,
                          int ncq, double tolerance,
                          spiro_batch_result *result);

/* Release all storage owned by result and clear it.  Safe for NULL and for an
   already-cleared result. */
void SpiroFreeBatchResult(spiro_batch_result *result);

#ifdef __cplusplus
}
#endif
#endif
