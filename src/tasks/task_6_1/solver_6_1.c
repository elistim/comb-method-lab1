#include <math.h>
#include <string.h>
#include "tasks/task_6_1/solver_6_1.h"

typedef struct {
    double left;
    double right;
    int valid;
} CandidateInterval;

/* Состояние комбинированного метода (Ньютон + хорда). */
typedef struct {
    double x_newton;
    double x_chord;
    double approx_hist[4];
    int hist_count;
} CombinedState;

/* Базовая валидация параметров поиска интервалов корней. */
static int is_separate_roots_args_valid(double (*func)(double),
                                        RootInterval *intervals,
                                        unsigned long *f_calls, double a,
                                        double b, double step, int max_count) {
    if (func == NULL || intervals == NULL || f_calls == NULL) {
        return 0;
    }
    if (a >= b || step <= 0.0 || max_count <= 0) {
        return 0;
    }
    return 1;
}

/* Выделяет кандидатный интервал: точный ноль или смена знака на шаге. */
static CandidateInterval detect_candidate_interval(double x_left,
                                                   double x_right,
                                                   double f_left,
                                                   double f_right, double tiny) {
    CandidateInterval c = {0.0, 0.0, 0};

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

/* Отсеивает повторные интервалы около одного и того же корня. */
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

/* Добавляет найденный интервал в массив результатов. */
static void append_interval(RootInterval *intervals, int *count,
                            CandidateInterval c) {
    intervals[*count].left = c.left;
    intervals[*count].right = c.right;
    (*count)++;
}

/* Инициализация статистики метода и значения корня "не определено". */
static void init_method_stats(MethodStats *stats, double *root) {
    memset(stats, 0, sizeof(*stats));
    stats->alpha = NAN;
    stats->last_dx = NAN;
    *root = NAN;
}

/* Нормализует границы: обеспечивает a <= b. */
static void normalize_bounds(double *a, double *b) {
    if (*a > *b) {
        double tmp = *a;
        *a = *b;
        *b = tmp;
    }
}

/* Обработка вырожденного интервала [a, a]. */
static int handle_degenerate_interval(double (*func)(double), double a,
                                      double eps2, double *root,
                                      MethodStats *stats) {
    double fa0 = func(a);
    stats->f_calls++;
    *root = a;
    if (fabs(fa0) <= eps2) {
        stats->converged = 1;
    }
    return 1;
}

/* Вычисляет значения на концах и проверяет базовые случаи останова. */
static int evaluate_bounds(double (*func)(double), double a, double b,
                           double *fa, double *fb, double *root,
                           MethodStats *stats) {
    *fa = func(a);
    *fb = func(b);
    stats->f_calls += 2;

    if (*fa == 0.0) {
        *root = a;
        stats->converged = 1;
        return 0;
    }
    if (*fb == 0.0) {
        *root = b;
        stats->converged = 1;
        return 0;
    }
    if ((*fa) * (*fb) > 0.0) {
        return 0;
    }
    return 1;
}

/* Выбирает, с какого конца делать шаг Ньютона, а с какого - шаг хордой. */
static void choose_start_points(double (*func_double_prime)(double), double a,
                                double b, double fa, double fb,
                                CombinedState *state, MethodStats *stats) {
    double fppa = func_double_prime(a);
    double fppb = func_double_prime(b);
    stats->fpp_calls += 2;

    if (fa * fppa > 0.0) {
        state->x_newton = a;
        state->x_chord = b;
    } else if (fb * fppb > 0.0) {
        state->x_newton = b;
        state->x_chord = a;
    } else {
        state->x_newton = a;
        state->x_chord = b;
    }
}

/* Обновляет историю последних 4 приближений. */
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

/* Выполняет одну итерацию комбинированного метода. */
static int run_combined_iteration(double (*func)(double),
                                  double (*func_prime)(double), double tiny,
                                  CombinedState *state, MethodStats *stats,
                                  double *x_curr, double *f_curr) {
    double fxn = func(state->x_newton);
    double fpxn = func_prime(state->x_newton);
    stats->f_calls++;
    stats->fp_calls++;

    if (fabs(fpxn) < tiny) {
        return 0;
    }
    double x_newton_next = state->x_newton - fxn / fpxn;

    double fxc = func(state->x_chord);
    double fxn_next = func(x_newton_next);
    stats->f_calls += 2;

    double denom = fxn_next - fxc;
    if (fabs(denom) < tiny) {
        return 0;
    }
    double x_chord_next =
        state->x_chord - fxc * (x_newton_next - state->x_chord) / denom;

    *x_curr = 0.5 * (x_newton_next + x_chord_next);
    *f_curr = func(*x_curr);
    stats->f_calls++;

    state->x_newton = x_newton_next;
    state->x_chord = x_chord_next;
    return 1;
}

/* Оценивает параметр alpha по истории сходимости, если данных достаточно. */
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
                   RootInterval *intervals, int max_count,
                   unsigned long *f_calls) {
    /* Допуск для сравнения с нулем и совпадения границ по double. */
    const double tiny = 1e-14;
    if (!is_separate_roots_args_valid(func, intervals, f_calls, a, b, step,
                                      max_count)) {
        return 0;
    }

    int count = 0;
    double x_left = a;
    double f_left = func(x_left);
    (*f_calls)++;

    /* Сканируем [a, b] шагом step и ищем интервалы с корнями. */
    while (x_left < b && count < max_count) {
        double x_right = x_left + step;
        if (x_right > b) {
            x_right = b;
        }

        double f_right = func(x_right);
        (*f_calls)++;

        /* Кандидатный интервал + фильтрация дублей соседних шагов. */
        CandidateInterval cand =
            detect_candidate_interval(x_left, x_right, f_left, f_right, tiny);
        if (cand.valid &&
            !is_duplicate_interval(intervals, count, cand, step, tiny) &&
            count < max_count) {
            append_interval(intervals, &count, cand);
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
    /* Ограничение числа итераций и защита от деления на почти ноль. */
    const int max_iter = 1000;
    const double tiny = 1e-15;
    CombinedState state = {0};
    double fa = 0.0;
    double fb = 0.0;

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

    if (!evaluate_bounds(func, a, b, &fa, &fb, root, stats)) {
        return;
    }

    choose_start_points(func_double_prime, a, b, fa, fb, &state, stats);
    push_history(&state, 0.5 * (state.x_newton + state.x_chord));

    /* Основной цикл: шаг Ньютона + шаг хордой, контроль критериев останова. */
    for (int iter = 1; iter <= max_iter; iter++) {
        double x_curr = 0.0;
        double f_curr = 0.0;
        double prev_approx = state.approx_hist[state.hist_count - 1];

        if (!run_combined_iteration(func, func_prime, tiny, &state, stats,
                                    &x_curr, &f_curr)) {
            break;
        }

        push_history(&state, x_curr);

        stats->iterations = iter;
        stats->last_dx = fabs(x_curr - prev_approx);

        /* Два критерия: сжатие отрезка и малость невязки. */
        if (fabs(state.x_newton - state.x_chord) <= eps1 && fabs(f_curr) <= eps2) {
            stats->converged = 1;
            break;
        }
    }

    *root = 0.5 * (state.x_newton + state.x_chord);
    estimate_alpha_if_possible(&state, tiny, stats);
}
