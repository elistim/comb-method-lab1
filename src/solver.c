#include <math.h>
#include <string.h>
#include "main.h"

int separate_roots(double (*func)(double), double a, double b, double step,
                   RootInterval *intervals, int max_count,
                   unsigned long *f_calls) {
    const double tiny = 1e-14;
    if (func == NULL || intervals == NULL || f_calls == NULL) {
        return 0;
    }
    if (a >= b || step <= 0.0 || max_count <= 0) {
        return 0;
    }

    int count = 0;
    double x_left = a;
    double f_left = func(x_left);
    (*f_calls)++;

    while (x_left < b && count < max_count) {
        double x_right = x_left + step;
        if (x_right > b) {
            x_right = b;
        }

        double f_right = func(x_right);
        (*f_calls)++;

        int should_add = 0;
        double cand_left = 0.0;
        double cand_right = 0.0;

        if (fabs(f_left) <= tiny) {
            cand_left = x_left;
            cand_right = x_left;
            should_add = 1;
        } else if (f_left * f_right < 0.0) {
            cand_left = x_left;
            cand_right = x_right;
            should_add = 1;
        } else if (fabs(f_right) <= tiny) {
            cand_left = x_right;
            cand_right = x_right;
            should_add = 1;
        }

        if (should_add) {
            int is_duplicate = 0;
            if (count > 0) {
                double prev_left = intervals[count - 1].left;
                double prev_right = intervals[count - 1].right;
                double prev_mid = 0.5 * (prev_left + prev_right);
                double cand_mid = 0.5 * (cand_left + cand_right);
                int same_endpoints =
                    fabs(prev_left - cand_left) <= tiny &&
                    fabs(prev_right - cand_right) <= tiny;
                int same_root_zone =
                    fabs(prev_mid - cand_mid) <= 0.5 * step &&
                    fabs(prev_right - prev_left) <= step &&
                    fabs(cand_right - cand_left) <= step;

                if (same_endpoints || same_root_zone) {
                    is_duplicate = 1;
                }
            }
            if (!is_duplicate && count < max_count) {
                intervals[count].left = cand_left;
                intervals[count].right = cand_right;
                count++;
            }
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
    double x_newton;
    double x_chord;
    double approx_hist[4] = {0.0, 0.0, 0.0, 0.0};
    int hist_count = 0;

    if (root == NULL || stats == NULL || func == NULL || func_prime == NULL ||
        func_double_prime == NULL) {
        return;
    }

    memset(stats, 0, sizeof(*stats));
    stats->alpha = NAN;
    stats->last_dx = NAN;
    *root = NAN;

    if (a > b) {
        double tmp = a;
        a = b;
        b = tmp;
    }
    if (eps1 <= 0.0 || eps2 <= 0.0) {
        return;
    }

    if (a == b) {
        double fa0 = func(a);
        stats->f_calls++;
        *root = a;
        if (fabs(fa0) <= eps2) {
            stats->converged = 1;
        }
        return;
    }

    double fa = func(a);
    double fb = func(b);
    stats->f_calls += 2;

    if (fa == 0.0) {
        *root = a;
        stats->converged = 1;
        return;
    }
    if (fb == 0.0) {
        *root = b;
        stats->converged = 1;
        return;
    }
    if (fa * fb > 0.0) {
        return;
    }

    double fppa = func_double_prime(a);
    double fppb = func_double_prime(b);
    stats->fpp_calls += 2;

    if (fa * fppa > 0.0) {
        x_newton = a;
        x_chord = b;
    } else if (fb * fppb > 0.0) {
        x_newton = b;
        x_chord = a;
    } else {
        x_newton = a;
        x_chord = b;
    }

    approx_hist[hist_count++] = 0.5 * (x_newton + x_chord);

    for (int iter = 1; iter <= max_iter; iter++) {
        double fxn = func(x_newton);
        double fpxn = func_prime(x_newton);
        stats->f_calls++;
        stats->fp_calls++;

        if (fabs(fpxn) < tiny) {
            break;
        }
        double x_newton_next = x_newton - fxn / fpxn;

        double fxc = func(x_chord);
        double fxn_next = func(x_newton_next);
        stats->f_calls += 2;

        double denom = fxn_next - fxc;
        if (fabs(denom) < tiny) {
            break;
        }
        double x_chord_next =
            x_chord - fxc * (x_newton_next - x_chord) / denom;

        double x_curr = 0.5 * (x_newton_next + x_chord_next);
        double f_curr = func(x_curr);
        stats->f_calls++;

        if (hist_count < 4) {
            approx_hist[hist_count++] = x_curr;
        } else {
            approx_hist[0] = approx_hist[1];
            approx_hist[1] = approx_hist[2];
            approx_hist[2] = approx_hist[3];
            approx_hist[3] = x_curr;
        }

        stats->iterations = iter;
        stats->last_dx = fabs(x_curr - approx_hist[hist_count - 2]);

        if (fabs(x_newton_next - x_chord_next) <= eps1 && fabs(f_curr) <= eps2) {
            x_newton = x_newton_next;
            x_chord = x_chord_next;
            stats->converged = 1;
            break;
        }

        x_newton = x_newton_next;
        x_chord = x_chord_next;
    }

    *root = 0.5 * (x_newton + x_chord);

    if (hist_count == 4) {
        double d0 = fabs(approx_hist[1] - approx_hist[0]);
        double d1 = fabs(approx_hist[2] - approx_hist[1]);
        double d2 = fabs(approx_hist[3] - approx_hist[2]);

        if (d0 > tiny && d1 > tiny && d2 > tiny && fabs(log(d1 / d0)) > tiny) {
            double k_est = log(d2 / d1) / log(d1 / d0);
            if (isfinite(k_est)) {
                stats->alpha = d2 / pow(d1, k_est);
            }
        }
    }
}
