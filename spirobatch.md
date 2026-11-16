# Batch conversion with error bounds

`spirobatch.h` adds a narrow batch entry point.  It converts several
independent contours in one call and returns flattened Bézier commands with an
absolute error bound on every geometric command.

## Interface

```c
typedef struct {
    const spiro_cp *points;
    int n;
    int is_closed;
} spiro_contour;

typedef struct {
    char ty;                    /* 'm', 'k', 'l', or 'c' */
    double x0, y0;
    double x1, y1;
    double x2, y2;
    double x3, y3;
    double err;                 /* absolute error upper bound */
} spiro_batch_segment;

typedef struct {
    spiro_batch_segment *segments;
    int segment_count;
    int segment_capacity;       /* internal */
    int *contour_offsets;       /* contour_count + 1 entries */
    double *contour_errors;     /* contour_count entries */
    int contour_count;
} spiro_batch_result;

int SpiroCPsToBezierTol(const spiro_cp *points, int n, int ncq,
                        int is_closed, double tolerance,
                        spiro_batch_result *result);

int SpiroCPsToBezierBatch(const spiro_contour *contours, int ncontours,
                          int ncq, double tolerance,
                          spiro_batch_result *result);

void SpiroFreeBatchResult(spiro_batch_result *result);
```

The commands in contour `c` are
`segments[contour_offsets[c] .. contour_offsets[c + 1])`.  `'m'` and `'k'`
carry no error; `'l'` and `'c'` have `err >= 0`.  `contour_errors[c]` is the
maximum geometric error in that contour.

```c
spiro_contour in[2] = {
    {points_a, count_a, 1},
    {points_b, count_b, 0}
};
spiro_batch_result out;

if (!SpiroCPsToBezierBatch(in, 2, 0, 0.01, &out))
    return failure;
/* consume out.segments / out.contour_offsets */
SpiroFreeBatchResult(&out);
```

`tolerance` is an absolute distance in input-coordinate units.  It must be
finite and greater than zero.  The single-contour and batch functions execute
the same conversion routine, so a contour in a batch is bit-for-bit identical
to converting it alone.

## Error bound

For every candidate cubic or line, the converter compares it with the exact
polynomial-spiral segment at 33 equally spaced parameters.  If
`d(t) = spiral(t) - bezier(t)`, the reported bound is

```text
max |d(t_i)| + h/2 * max |d'(t_i)|
             + h^2/8 * (max |spiral''| + max |bezier''|)
             + roundoff_floor
```

where `h = 1/32`.  The second-derivative terms bound the error between
samples; the last term is a conservative numerical-integration floor.  If the
bound is above `tolerance`, the segment is bisected.  Bisection is capped at
depth 14.  If the cap or the floating-point floor is reached, conversion fails
instead of returning a command whose claimed bound is false.

The error is a sup over the normalized Spiro parameter, not a Hausdorff
distance to a reparameterized curve.

## Boundary behaviour

| Input | Result |
| --- | --- |
| `ncontours == 0` | Success, empty result |
| `result == NULL` | Failure |
| `contours == NULL` with `ncontours > 0` | Failure |
| `points == NULL`, `n < 2`, or non-finite coordinate | Failure |
| `tolerance <= 0`, NaN, or infinity | Failure |
| Unsupported `ncq` bits | Failure |
| Failure in any contour | Atomic failure; result is empty |
| Very small tolerance below the roundoff floor | Failure, never an overstated bound |
| Open contour without `'{'` / `'}'` tags | Accepted; `is_closed` is authoritative |
| Closed contour without trailing `'z'` | Accepted; `is_closed` is authoritative |

Supported `ncq` bits are `SPIRO_RETRO_VER1` and
`SPIRO_INCLUDE_LAST_KNOT`.  `SPIRO_INTERNAL_BEZCTX`, `SPIRO_REVERSE_SRC`, and
the arc/quad output selectors are rejected because this interface always emits
cubic/line error-bounded output and never mutates the source.

## Memory

Contours are solved one at a time.  Solver memory is therefore proportional to
the largest contour, while the flattened command and offset arrays are
proportional to the total input point count for a fixed subdivision depth.
Peak storage is `O(total points + emitted commands)`, and emitted commands are
`O(total points)` because the depth cap is fixed.

## Consistency data

`tests/batch-test.c` uses three contours (16-point closed, 12-point open with
anchors/handles, 4-point closed circle; 32 points total), compares every batch
command with a standalone `SpiroCPsToBezierTol()` call, and checks every
`err <= tolerance`.

| tolerance | commands | geometric commands | max reported err |
| ---: | ---: | ---: | ---: |
| 1 | 75 | 45 | 0.750850150 |
| 0.1 | 109 | 79 | 0.0988569933 |
| 0.01 | 176 | 146 | 0.00965778374 |
| 0.001 | 324 | 294 | 0.00094564158 |
| 0.0001 | 721 | 691 | 9.95610205e-05 |

Batch and sequential results matched exactly for both `ncq = 0` and
`ncq = SPIRO_RETRO_VER1`.  Repeating the 16-point contour 16 times (256 input
points) produced 832 commands at `tolerance = 0.1`, exactly 16 times the
single-contour count.
