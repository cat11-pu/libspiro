#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "spiroentrypoints.h"

#define CASE_COUNT 200
#define MIN_CP 3
#define MAX_CP 64
#define MAX_CURVE 65536
#define HUGE_COORD 1.0e15
#define MAX_COUNTEREXAMPLES 16
#define TEST_PI 3.14159265358979323846

typedef enum {
    FEATURE_DUPLICATE,
    FEATURE_COLLINEAR,
    FEATURE_TINY_ANGLE,
    FEATURE_HUGE
} feature_kind;

typedef enum {
    OUTCOME_OK,
    OUTCOME_REJECTED,
    OUTCOME_NONFINITE,
    OUTCOME_OVERFLOW,
    OUTCOME_CRASH,
    OUTCOME_TIMEOUT,
    OUTCOME_NONDETERMINISTIC,
    OUTCOME_HARNESS_ERROR
} outcome_kind;

typedef struct {
    spiro_cp cp[MAX_CP];
    int n;
    feature_kind feature;
} test_case;

typedef struct {
    int32_t status;
    int32_t ret;
    int32_t count;
    int32_t bad_index;
    int32_t bad_field;
} run_header;

typedef struct {
    run_header header;
    curve_data *curves;
    int term_signal;
    int read_error;
} run_result;

typedef struct {
    outcome_kind kind;
} evaluation;

typedef struct {
    test_case tc;
    outcome_kind kind;
} counterexample;

static const char *feature_name(feature_kind feature) {
    switch (feature) {
    case FEATURE_DUPLICATE:
        return "duplicate";
    case FEATURE_COLLINEAR:
        return "collinear";
    case FEATURE_TINY_ANGLE:
        return "tiny-angle";
    case FEATURE_HUGE:
        return "huge";
    }
    return "unknown";
}

static const char *outcome_name(outcome_kind kind) {
    switch (kind) {
    case OUTCOME_OK:
        return "ok";
    case OUTCOME_REJECTED:
        return "rejected";
    case OUTCOME_NONFINITE:
        return "nonfinite";
    case OUTCOME_OVERFLOW:
        return "overflow";
    case OUTCOME_CRASH:
        return "crash";
    case OUTCOME_TIMEOUT:
        return "timeout";
    case OUTCOME_NONDETERMINISTIC:
        return "nondeterministic";
    case OUTCOME_HARNESS_ERROR:
        return "harness-error";
    }
    return "unknown";
}

static uint64_t rng_next(uint64_t *state) {
    uint64_t z;

    *state += UINT64_C(0x9e3779b97f4a7c15);
    z = *state;
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static double rng_unit(uint64_t *state) {
    return (double)(rng_next(state) >> 11) * (1.0 / 9007199254740992.0);
}

static void make_glyph_case(int id, test_case *tc) {
    uint64_t state = UINT64_C(0x6a09e667f3bcc909) ^ (uint64_t)(id + 1);
    double cx = 500.0 + 180.0 * rng_unit(&state);
    double cy = 650.0 + 220.0 * rng_unit(&state);
    double rx = 220.0 + 160.0 * rng_unit(&state);
    double ry = 300.0 + 240.0 * rng_unit(&state);
    double skew = 0.12 + 0.20 * rng_unit(&state);
    int n = 18 + (id % 7);
    int i;
    int k;

    memset(tc, 0, sizeof(*tc));
    tc->n = n;
    tc->feature = (feature_kind)(id % 4);

    for (i = 0; i < n; i++) {
        double theta = (2.0 * TEST_PI * (double)i) / (double)n;
        double c = cos(theta);
        double s = sin(theta);
        double ex = copysign(pow(fabs(c), 0.62), c);
        double ey = copysign(pow(fabs(s), 0.72), s);
        double jitter = 1.0 + 0.035 * (2.0 * rng_unit(&state) - 1.0);

        tc->cp[i].x = cx + rx * ex * jitter + skew * ry * ey;
        tc->cp[i].y = cy + ry * ey * jitter - skew * rx * ex;
        if (i == 0 || i == n - 1 || i % 5 == 0)
            tc->cp[i].ty = SPIRO_CORNER;
        else if (i % 3 == 0)
            tc->cp[i].ty = SPIRO_G4;
        else
            tc->cp[i].ty = SPIRO_G2;
    }
    tc->cp[0].ty = SPIRO_CORNER;
    tc->cp[n - 1].ty = SPIRO_CORNER;

    switch (tc->feature) {
    case FEATURE_DUPLICATE:
        k = 2 + (id % (n - 5));
        tc->cp[k + 1] = tc->cp[k];
        tc->cp[k + 1].ty = SPIRO_CORNER;
        break;
    case FEATURE_COLLINEAR:
        k = 2 + ((id * 3) % (n - 5));
        tc->cp[k + 1].x = 0.5 * (tc->cp[k].x + tc->cp[k + 2].x);
        tc->cp[k + 1].y = 0.5 * (tc->cp[k].y + tc->cp[k + 2].y);
        tc->cp[k + 1].ty = SPIRO_G2;
        break;
    case FEATURE_TINY_ANGLE: {
        double dx;
        double dy;
        double length;
        double nx;
        double ny;
        double epsilon;

        k = 2 + ((id * 5) % (n - 5));
        dx = tc->cp[k + 2].x - tc->cp[k].x;
        dy = tc->cp[k + 2].y - tc->cp[k].y;
        length = hypot(dx, dy);
        if (length == 0.0)
            length = 1.0;
        nx = -dy / length;
        ny = dx / length;
        epsilon = 1.0e-11 * (1.0 + (double)(id % 9));
        tc->cp[k + 1].x = 0.5 * (tc->cp[k].x + tc->cp[k + 2].x) + nx * epsilon;
        tc->cp[k + 1].y = 0.5 * (tc->cp[k].y + tc->cp[k + 2].y) + ny * epsilon;
        tc->cp[k + 1].ty = SPIRO_G2;
        break;
    }
    case FEATURE_HUGE:
        for (i = 3; i < n - 3; i += 2) {
            double sign = ((i + id) & 1) ? 1.0 : -1.0;
            tc->cp[i].x = sign * HUGE_COORD + tc->cp[i].x;
            tc->cp[i].y = -sign * HUGE_COORD + tc->cp[i].y;
        }
        break;
    }

}

static int write_all(int fd, const void *buffer, size_t length) {
    const unsigned char *p = (const unsigned char *)buffer;

    while (length > 0) {
        ssize_t written = write(fd, p, length);
        if (written < 0) {
            if (errno == EINTR)
                continue;
            return 0;
        }
        p += (size_t)written;
        length -= (size_t)written;
    }
    return 1;
}

static int write_i32(int fd, int32_t value) {
    return write_all(fd, &value, sizeof(value));
}

static int write_u64(int fd, uint64_t value) {
    return write_all(fd, &value, sizeof(value));
}

static uint64_t double_bits(double value) {
    uint64_t bits;

    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

static void child_evaluate(const test_case *tc, int fd) {
    spiro_cp cp[MAX_CP];
    ls_bezctx *bc;
    run_header header;
    int ret;
    int i;
    int j;
    int count;
    int status = OUTCOME_OK;
    int bad_index = -1;
    int bad_field = -1;
    static const int ncq = SPIRO_CUBIC_TO_BEZIER | SPIRO_INTERNAL_BEZCTX;

    memset(&header, 0, sizeof(header));
    memcpy(cp, tc->cp, (size_t)tc->n * sizeof(cp[0]));
    bc = new_ls_bezctx(MAX_CURVE, ncq);
    if (bc == NULL) {
        header.status = OUTCOME_HARNESS_ERROR;
        header.ret = 0;
        header.count = -1;
        write_i32(fd, header.status);
        write_i32(fd, header.ret);
        write_i32(fd, header.count);
        write_i32(fd, header.bad_index);
        write_i32(fd, header.bad_field);
        _exit(0);
    }

    ret = SpiroCPsToBezier2(cp, tc->n, ncq, 1, (bezctx *)bc);
    count = bc->l;
    if (ret == 0)
        status = OUTCOME_REJECTED;
    else if (count < 0)
        status = OUTCOME_OVERFLOW;
    else {
        for (i = 0; i < count && status == OUTCOME_OK; i++) {
            const double values[8] = {
                bc->cd[i].x0, bc->cd[i].y0,
                bc->cd[i].x1, bc->cd[i].y1,
                bc->cd[i].x2, bc->cd[i].y2,
                bc->cd[i].x3, bc->cd[i].y3
            };
            for (j = 0; j < 8; j++) {
                if (!isfinite(values[j])) {
                    status = OUTCOME_NONFINITE;
                    bad_index = i;
                    bad_field = j;
                    break;
                }
            }
        }
    }

    header.status = status;
    header.ret = ret;
    header.count = count;
    header.bad_index = bad_index;
    header.bad_field = bad_field;
    if (!write_i32(fd, header.status) ||
        !write_i32(fd, header.ret) ||
        !write_i32(fd, header.count) ||
        !write_i32(fd, header.bad_index) ||
        !write_i32(fd, header.bad_field))
        _exit(2);

    if (count > 0 && count <= MAX_CURVE) {
        for (i = 0; i < count; i++) {
            int32_t ty = (int32_t)bc->cd[i].ty;
            const double values[8] = {
                bc->cd[i].x0, bc->cd[i].y0,
                bc->cd[i].x1, bc->cd[i].y1,
                bc->cd[i].x2, bc->cd[i].y2,
                bc->cd[i].x3, bc->cd[i].y3
            };
            if (!write_i32(fd, ty))
                _exit(2);
            for (j = 0; j < 8; j++)
                if (!write_u64(fd, double_bits(values[j])))
                    _exit(2);
        }
    }
    free_ls_bezctx(bc);
    _exit(0);
}

static int grow_buffer(unsigned char **buffer, size_t *capacity, size_t needed) {
    unsigned char *next;
    size_t next_capacity;

    if (needed <= *capacity)
        return 1;
    next_capacity = *capacity ? *capacity : 4096;
    while (next_capacity < needed)
        next_capacity *= 2;
    next = (unsigned char *)realloc(*buffer, next_capacity);
    if (next == NULL)
        return 0;
    *buffer = next;
    *capacity = next_capacity;
    return 1;
}

static run_result run_once(const test_case *tc) {
    run_result result;
    unsigned char *buffer = NULL;
    size_t length = 0;
    size_t capacity = 0;
    int fds[2];
    pid_t pid;
    int wait_status = 0;
    size_t offset = 0;
    int i;

    memset(&result, 0, sizeof(result));
    result.header.status = OUTCOME_HARNESS_ERROR;
    result.header.count = -1;
    result.term_signal = 0;
    result.read_error = 0;

    if (pipe(fds) != 0) {
        result.read_error = 1;
        return result;
    }
    pid = fork();
    if (pid < 0) {
        close(fds[0]);
        close(fds[1]);
        result.read_error = 1;
        return result;
    }
    if (pid == 0) {
        close(fds[0]);
        alarm(15);
        child_evaluate(tc, fds[1]);
        _exit(2);
    }

    close(fds[1]);
    for (;;) {
        ssize_t got;

        if (!grow_buffer(&buffer, &capacity, length + 4096)) {
            result.read_error = 1;
            break;
        }
        got = read(fds[0], buffer + length, capacity - length);
        if (got < 0) {
            if (errno == EINTR)
                continue;
            result.read_error = 1;
            break;
        }
        if (got == 0)
            break;
        length += (size_t)got;
    }
    close(fds[0]);
    while (waitpid(pid, &wait_status, 0) < 0) {
        if (errno != EINTR) {
            result.read_error = 1;
            break;
        }
    }

    if (WIFSIGNALED(wait_status)) {
        result.term_signal = WTERMSIG(wait_status);
        if (result.term_signal == SIGALRM)
            result.header.status = OUTCOME_TIMEOUT;
        else
            result.header.status = OUTCOME_CRASH;
        free(buffer);
        return result;
    }
    if (!WIFEXITED(wait_status) || WEXITSTATUS(wait_status) != 0 ||
        length < 5 * sizeof(int32_t)) {
        result.read_error = 1;
        free(buffer);
        return result;
    }

    memcpy(&result.header.status, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);
    memcpy(&result.header.ret, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);
    memcpy(&result.header.count, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);
    memcpy(&result.header.bad_index, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);
    memcpy(&result.header.bad_field, buffer + offset, sizeof(int32_t));
    offset += sizeof(int32_t);

    if (result.header.count > 0 && result.header.count <= MAX_CURVE) {
        size_t row_size = sizeof(int32_t) + 8 * sizeof(uint64_t);
        size_t data_size = (size_t)result.header.count * row_size;

        if (length < offset + data_size) {
            result.read_error = 1;
            free(buffer);
            return result;
        }
        result.curves = (curve_data *)calloc((size_t)result.header.count,
                                             sizeof(curve_data));
        if (result.curves == NULL) {
            result.read_error = 1;
            free(buffer);
            return result;
        }
        for (i = 0; i < result.header.count; i++) {
            int32_t ty;
            double values[8];
            int j;

            memcpy(&ty, buffer + offset, sizeof(ty));
            offset += sizeof(ty);
            for (j = 0; j < 8; j++) {
                uint64_t bits;
                memcpy(&bits, buffer + offset, sizeof(bits));
                offset += sizeof(bits);
                memcpy(&values[j], &bits, sizeof(double));
            }
            result.curves[i].ty = (char)ty;
            result.curves[i].x0 = values[0];
            result.curves[i].y0 = values[1];
            result.curves[i].x1 = values[2];
            result.curves[i].y1 = values[3];
            result.curves[i].x2 = values[4];
            result.curves[i].y2 = values[5];
            result.curves[i].x3 = values[6];
            result.curves[i].y3 = values[7];
        }
    }
    free(buffer);
    return result;
}

static int same_run(const run_result *a, const run_result *b) {
    int i;

    if (a->term_signal != b->term_signal || a->read_error != b->read_error)
        return 0;
    if (a->header.status != b->header.status ||
        a->header.ret != b->header.ret ||
        a->header.count != b->header.count ||
        a->header.bad_index != b->header.bad_index ||
        a->header.bad_field != b->header.bad_field)
        return 0;
    if (a->header.count > 0 && a->header.count <= MAX_CURVE) {
        for (i = 0; i < a->header.count; i++) {
            const double av[8] = {
                a->curves[i].x0, a->curves[i].y0,
                a->curves[i].x1, a->curves[i].y1,
                a->curves[i].x2, a->curves[i].y2,
                a->curves[i].x3, a->curves[i].y3
            };
            const double bv[8] = {
                b->curves[i].x0, b->curves[i].y0,
                b->curves[i].x1, b->curves[i].y1,
                b->curves[i].x2, b->curves[i].y2,
                b->curves[i].x3, b->curves[i].y3
            };
            int j;

            if (a->curves[i].ty != b->curves[i].ty)
                return 0;
            for (j = 0; j < 8; j++)
                if (double_bits(av[j]) != double_bits(bv[j]))
                    return 0;
        }
    }
    return 1;
}

static evaluation evaluate(const test_case *tc) {
    evaluation eval;
    run_result first = run_once(tc);
    run_result second = run_once(tc);
    int deterministic;

    memset(&eval, 0, sizeof(eval));
    deterministic = same_run(&first, &second);

    if (first.term_signal == SIGALRM || second.term_signal == SIGALRM) {
        eval.kind = OUTCOME_TIMEOUT;
    } else if (first.term_signal != 0 || second.term_signal != 0) {
        eval.kind = OUTCOME_CRASH;
    } else if (first.read_error || second.read_error) {
        eval.kind = OUTCOME_HARNESS_ERROR;
    } else if (!deterministic) {
        eval.kind = OUTCOME_NONDETERMINISTIC;
    } else if (first.header.status != OUTCOME_OK) {
        eval.kind = (outcome_kind)first.header.status;
    } else {
        eval.kind = OUTCOME_OK;
    }

    free(first.curves);
    free(second.curves);
    return eval;
}

static int minimize_case(const test_case *original, test_case *minimal,
                         outcome_kind target) {
    test_case current = *original;
    int granularity = 2;
    int changed;

    for (;;) {
        int chunk;
        int start;

        changed = 0;
        if (current.n <= MIN_CP)
            break;
        chunk = (current.n + granularity - 1) / granularity;
        if (chunk < 1)
            chunk = 1;
        for (start = 0; start < current.n; start += chunk) {
            test_case candidate;
            int end = start + chunk;
            int keep;

            if (end > current.n)
                end = current.n;
            if (current.n - (end - start) < MIN_CP)
                continue;
            memset(&candidate, 0, sizeof(candidate));
            candidate.feature = current.feature;
            candidate.n = 0;
            for (keep = 0; keep < current.n; keep++) {
                if (keep >= start && keep < end)
                    continue;
                candidate.cp[candidate.n++] = current.cp[keep];
            }
            if (candidate.n < MIN_CP)
                continue;
            if (evaluate(&candidate).kind == target) {
                current = candidate;
                changed = 1;
                if (granularity > 2)
                    granularity--;
                break;
            }
        }
        if (!changed) {
            if (granularity >= current.n)
                break;
            granularity *= 2;
            if (granularity > current.n)
                granularity = current.n;
        }
    }

    for (;;) {
        int index;
        int removed = 0;

        if (current.n <= MIN_CP)
            break;
        for (index = 0; index < current.n; index++) {
            test_case candidate;
            int i;

            memset(&candidate, 0, sizeof(candidate));
            candidate.feature = current.feature;
            for (i = 0; i < current.n; i++)
                if (i != index)
                    candidate.cp[candidate.n++] = current.cp[i];
            if (evaluate(&candidate).kind == target) {
                current = candidate;
                removed = 1;
                break;
            }
        }
        if (!removed)
            break;
    }
    *minimal = current;
    return 1;
}

static void print_case(const test_case *tc) {
    int i;

    printf("spiro_cp minimal[] = {\n");
    for (i = 0; i < tc->n; i++) {
        printf("    {%.17g, %.17g, '%c'}%s\n",
               tc->cp[i].x, tc->cp[i].y, tc->cp[i].ty,
               i + 1 == tc->n ? "" : ",");
    }
    printf("};\n");
    printf("n=%d closed=1 ncq=0x%x\n\n", tc->n,
           SPIRO_CUBIC_TO_BEZIER | SPIRO_INTERNAL_BEZCTX);
}

static void usage(const char *argv0) {
    fprintf(stderr,
            "usage: %s [--count N] [--no-minimize] [--all] [--require-success]\n",
            argv0);
}

int main(int argc, char **argv) {
    int count = CASE_COUNT;
    int no_minimize = 0;
    int show_all = 0;
    int require_success = 0;
    int i;
    int failures = 0;
    int crashes = 0;
    int nonfinite = 0;
    int nondeterministic = 0;
    int rejected = 0;
    int overflow = 0;
    int timeouts = 0;
    int harness_errors = 0;
    int hard_failures = 0;
    int feature_total[4] = {0, 0, 0, 0};
    int feature_fail[4] = {0, 0, 0, 0};
    counterexample unique[64];
    int unique_count = 0;

    for (i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--count") == 0 && i + 1 < argc) {
            count = atoi(argv[++i]);
            if (count < 0 || count > 100000) {
                usage(argv[0]);
                return 2;
            }
        } else if (strcmp(argv[i], "--no-minimize") == 0) {
            no_minimize = 1;
        } else if (strcmp(argv[i], "--all") == 0) {
            show_all = 1;
        } else if (strcmp(argv[i], "--require-success") == 0) {
            require_success = 1;
        } else {
            usage(argv[0]);
            return 2;
        }
    }

    printf("libspiro property test: cases=%d mode=closed/cubic seed=fixed\n", count);
    for (i = 0; i < count; i++) {
        test_case tc;
        evaluation eval;

        make_glyph_case(i, &tc);
        feature_total[(int)tc.feature]++;
        eval = evaluate(&tc);
        if (eval.kind == OUTCOME_OK)
            continue;

        failures++;
        feature_fail[(int)tc.feature]++;
        switch (eval.kind) {
        case OUTCOME_CRASH:
            crashes++;
            break;
        case OUTCOME_TIMEOUT:
            timeouts++;
            break;
        case OUTCOME_NONFINITE:
            nonfinite++;
            break;
        case OUTCOME_NONDETERMINISTIC:
            nondeterministic++;
            break;
        case OUTCOME_REJECTED:
            rejected++;
            break;
        case OUTCOME_OVERFLOW:
            overflow++;
            break;
        default:
            harness_errors++;
            break;
        }
        if (eval.kind == OUTCOME_CRASH || eval.kind == OUTCOME_TIMEOUT ||
            eval.kind == OUTCOME_NONFINITE ||
            eval.kind == OUTCOME_OVERFLOW ||
            eval.kind == OUTCOME_NONDETERMINISTIC ||
            eval.kind == OUTCOME_HARNESS_ERROR)
            hard_failures++;
        if (!no_minimize &&
            (unique_count < (int)(sizeof(unique) / sizeof(unique[0])))) {
            test_case minimal;
            int found = -1;
            int j;

            minimize_case(&tc, &minimal, eval.kind);
            for (j = 0; j < unique_count; j++) {
                if (unique[j].kind == eval.kind &&
                    unique[j].tc.feature == minimal.feature) {
                    found = j;
                    break;
                }
            }
            if (found >= 0) {
                if (minimal.n < unique[found].tc.n)
                    unique[found].tc = minimal;
            } else {
                unique[unique_count].tc = minimal;
                unique[unique_count].kind = eval.kind;
                unique_count++;
            }
        } else if (no_minimize &&
                   unique_count < (int)(sizeof(unique) / sizeof(unique[0]))) {
            unique[unique_count].tc = tc;
            unique[unique_count].kind = eval.kind;
            unique_count++;
        }
    }

    printf("\nsummary\n");
    printf("  cases=%d runs=%d\n", count, count * 2);
    printf("  accepted=%d rejected=%d overflow=%d\n",
           count - failures, rejected, overflow);
    printf("  crashes=%d nonfinite=%d nondeterministic=%d timeout=%d "
           "harness_errors=%d\n",
           crashes, nonfinite, nondeterministic, timeouts, harness_errors);
    printf("  hard_property_violations=%d\n", hard_failures);
    printf("  feature totals/failures: duplicate=%d/%d collinear=%d/%d "
           "tiny-angle=%d/%d huge=%d/%d\n",
           feature_total[0], feature_fail[0],
           feature_total[1], feature_fail[1],
           feature_total[2], feature_fail[2],
           feature_total[3], feature_fail[3]);

    if (!no_minimize) {
        int shown = unique_count;
        int limit = show_all ? unique_count : MAX_COUNTEREXAMPLES;

        if (shown > limit)
            shown = limit;
        printf("\nrejection counterexamples: unique_minimal=%d shown=%d\n",
               unique_count, shown);
        for (i = 0; i < shown; i++) {
            printf("\n[%d] outcome=%s feature=%s\n", i + 1,
                   outcome_name(unique[i].kind),
                   feature_name(unique[i].tc.feature));
            print_case(&unique[i].tc);
        }
        if (unique_count > shown)
            printf("... %d more omitted; use --all to print all\n",
                   unique_count - shown);
    }

    if (hard_failures != 0)
        return 1;
    if (require_success && failures != 0)
        return 1;
    return 0;
}
