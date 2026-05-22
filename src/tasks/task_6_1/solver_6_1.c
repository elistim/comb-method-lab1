#include <math.h>
#include <string.h>
#include "tasks/task_6_1/solver_6_1.h"

typedef struct {
    double left;
    double right;
    int valid;
} CandidateInterval;

typedef struct {
    double a;
    double b;
    double fa;
    double fb;
    double approx_hist[4];
    int hist_count;
} CombinedState;

static int is_separate_roots_args_valid(double (*func)(double),
                                        RootInterval *intervals,
                                        unsigned long *f_calls, double a,
                                        double b, double step, int capacity) {
    return func != NULL && intervals != NULL && f_calls != NULL && a < b &&
           step > 0.0 && capacity > 0;
}

static CandidateInterval detect_candidate_interval(double x_left,
                                                   double x_right,
                                                   double f_left,
                                                   double f_right,
                                                   double tiny) {
    CandidateInterval c = {0.0, 0.0, 0};

    if (!isfinite(f_left) || !isfinite(f_right)) {
        return c;
    }

    if (fabs(f_left) <= tiny) {
        c.left = x_left;
        c.right = x_left;
        c.valid = 1;
    } else if (f_left * f_right < 0.0) {
        c.left = x_left;
        c.right = x_right;
        c.valid = 1;
    } else if (fabs(f_right) <= tiny) {
        c.left = x_right;
        c.right = x_right;
        c.valid = 1;
    }

    return c;
}

static int is_duplicate_interval(const RootInterval *intervals, int count,
                                 CandidateInterval c, double step,
                                 double tiny) {
    if (count <= 0) {
        return 0;
    }

    double prev_left = intervals[count - 1].left;
    double prev_right = intervals[count - 1].right;
    double prev_mid = 0.5 * (prev_left + prev_right);
    double cand_mid = 0.5 * (c.left + c.right);
    int same_endpoints =
        fabs(prev_left - c.left) <= tiny && fabs(prev_right - c.right) <= tiny;
    int same_root_zone = fabs(prev_mid - cand_mid) <= 0.5 * step &&
                         fabs(prev_right - prev_left) <= step &&
                         fabs(c.right - c.left) <= step;

    return same_endpoints || same_root_zone;
}

static void init_method_stats(MethodStats *stats, double *root) {
    memset(stats, 0, sizeof(*stats));
    stats->alpha = NAN;
    stats->last_dx = NAN;
    stats->argument_error = NAN;
    stats->residual_abs = NAN;
    *root = NAN;
}

static void normalize_bounds(double *a, double *b) {
    if (*a > *b) {
        double tmp = *a;
        *a = *b;
        *b = tmp;
    }
}

static int handle_degenerate_interval(double (*func)(double), double a,
                                      double eps2, double *root,
                                      MethodStats *stats) {
    double fa0 = func(a);
    stats->f_calls++;
    *root = a;
    stats->argument_error = 0.0;
    stats->residual_abs = fabs(fa0);
    if (fabs(fa0) <= eps2) {
        stats->converged = 1;
    }
    return 1;
}

static int evaluate_bounds(double (*func)(double), double a, double b,
                           CombinedState *state, double *root,
                           MethodStats *stats, double tiny) {
    state->a = a;
    state->b = b;
    state->fa = func(a);
    state->fb = func(b);
    stats->f_calls += 2;

    if (!isfinite(state->fa) || !isfinite(state->fb)) {
        return 0;
    }
    if (fabs(state->fa) <= tiny) {
        *root = a;
        stats->argument_error = 0.0;
        stats->residual_abs = fabs(state->fa);
        stats->converged = 1;
        return 0;
    }
    if (fabs(state->fb) <= tiny) {
        *root = b;
        stats->argument_error = 0.0;
        stats->residual_abs = fabs(state->fb);
        stats->converged = 1;
        return 0;
    }
    if (state->fa * state->fb > 0.0) {
        return 0;
    }
    return 1;
}

static void push_history(CombinedState *state, double x_curr) {
    if (state->hist_count < 4) {
        state->approx_hist[state->hist_count++] = x_curr;
        return;
    }

    state->approx_hist[0] = state->approx_hist[1];
    state->approx_hist[1] = state->approx_hist[2];
    state->approx_hist[2] = state->approx_hist[3];
    state->approx_hist[3] = x_curr;
}

static int choose_newton_side(double (*func_double_prime)(double),
                              const CombinedState *state, int *from_left,
                              MethodStats *stats, double tiny) {
    double fppa = func_double_prime(state->a);
    double fppb = func_double_prime(state->b);
    stats->fpp_calls += 2;

    if (!isfinite(fppa) || !isfinite(fppb)) {
        return 0;
    }
    if (state->fa * fppa > 0.0) {
        *from_left = 1;
        return 1;
    }
    if (state->fb * fppb > 0.0) {
        *from_left = 0;
        return 1;
    }

    *from_left = fabs(state->fa * fppa) >= fabs(state->fb * fppb);
    return fabs(fppa) > tiny || fabs(fppb) > tiny;
}

static int run_combined_iteration(double (*func)(double),
                                  double (*func_prime)(double),
                                  double (*func_double_prime)(double),
                                  double tiny, CombinedState *state,
                                  MethodStats *stats) {
    int newton_from_left = 1;
    double x_newton = 0.0;
    double x_chord = 0.0;
    double f_newton = 0.0;
    double f_chord = 0.0;
    double x0 = 0.0;
    double f0 = 0.0;
    double fp0 = 0.0;
    double denom = state->fb - state->fa;

    if (fabs(denom) < tiny ||
        !choose_newton_side(func_double_prime, state, &newton_from_left, stats,
                            tiny)) {
        return 0;
    }

    if (newton_from_left) {
        x0 = state->a;
        f0 = state->fa;
    } else {
        x0 = state->b;
        f0 = state->fb;
    }

    fp0 = func_prime(x0);
    stats->fp_calls++;
    if (!isfinite(fp0) || fabs(fp0) < tiny) {
        return 0;
    }

    x_newton = x0 - f0 / fp0;
    x_chord = state->a - state->fa * (state->b - state->a) / denom;

    if (x_newton < state->a || x_newton > state->b ||
        x_chord < state->a || x_chord > state->b) {
        return 0;
    }

    f_newton = func(x_newton);
    f_chord = func(x_chord);
    stats->f_calls += 2;
    if (!isfinite(f_newton) || !isfinite(f_chord)) {
        return 0;
    }

    if (x_newton <= x_chord) {
        state->a = x_newton;
        state->fa = f_newton;
        state->b = x_chord;
        state->fb = f_chord;
    } else {
        state->a = x_chord;
        state->fa = f_chord;
        state->b = x_newton;
        state->fb = f_newton;
    }

    return 1;
}

static void estimate_alpha_if_possible(const CombinedState *state, double tiny,
                                       MethodStats *stats) {
    if (state->hist_count != 4) {
        return;
    }

    double d0 = fabs(state->approx_hist[1] - state->approx_hist[0]);
    double d1 = fabs(state->approx_hist[2] - state->approx_hist[1]);
    double d2 = fabs(state->approx_hist[3] - state->approx_hist[2]);

    if (d0 > tiny && d1 > tiny && d2 > tiny && fabs(log(d1 / d0)) > tiny) {
        double k_est = log(d2 / d1) / log(d1 / d0);
        if (isfinite(k_est)) {
            stats->alpha = d2 / pow(d1, k_est);
        }
    }
}

int separate_roots(double (*func)(double), double a, double b, double step,
                   RootInterval *intervals, int capacity,
                   unsigned long *f_calls) {
    const double tiny = 1e-14;
    int count = 0;
    double x_left = a;
    double f_left = 0.0;

    if (!is_separate_roots_args_valid(func, intervals, f_calls, a, b, step,
                                      capacity)) {
        return 0;
    }

    f_left = func(x_left);
    (*f_calls)++;

    while (x_left < b) {
        double x_right = x_left + step;
        double f_right = 0.0;
        CandidateInterval cand;

        if (x_right > b) {
            x_right = b;
        }

        f_right = func(x_right);
        (*f_calls)++;

        cand =
            detect_candidate_interval(x_left, x_right, f_left, f_right, tiny);
        if (cand.valid &&
            !is_duplicate_interval(intervals, count, cand, step, tiny)) {
            if (count >= capacity) {
                break;
            }
            intervals[count].left = cand.left;
            intervals[count].right = cand.right;
            count++;
        }

        x_left = x_right;
        f_left = f_right;
    }

    return count;
}

void combined_method(double (*func)(double), double (*func_prime)(double),
                     double (*func_double_prime)(double), double a, double b,
                     double eps1, double eps2, double *root,
                     MethodStats *stats) {
    const int max_iter = 1000;
    const double tiny = 1e-15;
    CombinedState state = {0};

    if (root == NULL || stats == NULL || func == NULL || func_prime == NULL ||
        func_double_prime == NULL) {
        return;
    }

    init_method_stats(stats, root);
    normalize_bounds(&a, &b);

    if (eps1 <= 0.0 || eps2 <= 0.0) {
        return;
    }

    if (a == b && handle_degenerate_interval(func, a, eps2, root, stats)) {
        return;
    }

    if (!evaluate_bounds(func, a, b, &state, root, stats, tiny)) {
        return;
    }

    push_history(&state, 0.5 * (state.a + state.b));

    for (int iter = 1; iter <= max_iter; iter++) {
        double prev_approx = state.approx_hist[state.hist_count - 1];
        double x_curr = 0.0;
        double f_curr = 0.0;

        if (!run_combined_iteration(func, func_prime, func_double_prime, tiny,
                                    &state, stats)) {
            break;
        }

        x_curr = 0.5 * (state.a + state.b);
        f_curr = fabs(state.fa) < fabs(state.fb) ? state.fa : state.fb;
        push_history(&state, x_curr);

        stats->iterations = iter;
        stats->last_dx = fabs(x_curr - prev_approx);
        stats->argument_error = 0.5 * fabs(state.b - state.a);
        stats->residual_abs = fabs(f_curr);

        if (stats->argument_error <= eps1 && fabs(f_curr) <= eps2) {
            stats->converged = 1;
            break;
        }
    }

    *root = 0.5 * (state.a + state.b);
    estimate_alpha_if_possible(&state, tiny, stats);
}
