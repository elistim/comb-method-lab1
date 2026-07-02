#include <math.h>
#include <stdio.h>
#include "scientific_format.h"
#include "tasks/task_6_2/task_6_2_1.h"

#define N 5

static void print_system(const double a[N][N], const double b[N]) {
    int i = 0;
    int j = 0;

    printf("Система A*x=b (вариант 9):\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%10.2f ", a[i][j]);
        }
        printf("| %8.2f\n", b[i]);
    }
    printf("\n");
}

static int is_strict_diagonal_dominance(const double a[N][N]) {
    int i = 0;
    int j = 0;

    for (i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (j = 0; j < N; j++) {
            if (j != i) {
                row_sum += fabs(a[i][j]);
            }
        }
        if (fabs(a[i][i]) <= row_sum) {
            return 0;
        }
    }

    return 1;
}

static void calc_residual(const double a[N][N], const double b[N],
                          const double x[N], double r[N],
                          double *residual_inf_norm) {
    int i = 0;
    int j = 0;

    *residual_inf_norm = 0.0;
    for (i = 0; i < N; i++) {
        double ax_i = 0.0;
        for (j = 0; j < N; j++) {
            ax_i += a[i][j] * x[j];
        }
        r[i] = ax_i - b[i];
        if (fabs(r[i]) > *residual_inf_norm) {
            *residual_inf_norm = fabs(r[i]);
        }
    }
}

static double calc_alpha_norm(const double a[N][N], double alpha_out[N][N],
                                double row_sums_out[N], int *max_row_out) {
    double q = 0.0;
    int max_row = 0;
    int i = 0;
    int j = 0;

    for (i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (j = 0; j < N; j++) {
            double alpha_ij = (i == j) ? 0.0 : -a[i][j] / a[i][i];
            row_sum += fabs(alpha_ij);
            if (alpha_out != NULL) {
                alpha_out[i][j] = alpha_ij;
            }
        }
        if (row_sums_out != NULL) {
            row_sums_out[i] = row_sum;
        }
        if (row_sum > q) {
            q = row_sum;
            max_row = i;
        }
    }

    if (max_row_out != NULL) {
        *max_row_out = max_row;
    }

    return q;
}

static void print_alpha_details(const double a[N][N], const double b[N],
                                const double alpha[N][N],
                                const double row_sums[N],
                                int max_row, double q) {
    int i = 0;
    int j = 0;

    printf("\nРасчёт коэффициентов приведённой системы x = beta + alpha*x:\n");
    printf("Коэффициенты beta_i = b_i / a_ii:\n");
    for (i = 0; i < N; i++) {
        printf("beta_%d = %.2f / %.2f = %.12f\n",
               i + 1, b[i], a[i][i], b[i] / a[i][i]);
    }

    printf("\nКоэффициенты alpha_ij = -a_ij / a_ii при i != j; alpha_ii = 0:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (i == j) {
                printf("alpha_%d%d = 0\n", i + 1, j + 1);
            } else {
                printf("alpha_%d%d = -(%.2f) / %.2f = %.12f\n",
                       i + 1, j + 1, a[i][j], a[i][i], alpha[i][j]);
            }
        }
    }

    printf("\nПолученная приведённая система:\n");
    for (i = 0; i < N; i++) {
        printf("x_%d = %.12f", i + 1, b[i] / a[i][i]);
        for (j = 0; j < N; j++) {
            if (j != i) {
                printf(" + (% .6f)*x_%d", alpha[i][j], j + 1);
            }
        }
        printf("\n");
    }

    printf("\nМатрица коэффициентов alpha:\n");
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("% .9f%s", alpha[i][j], (j + 1 < N) ? "  " : "");
        }
        printf("\n");
    }

    printf("\nСуммы модулей коэффициентов по строкам alpha:\n");
    for (i = 0; i < N; i++) {
        printf("Строка %d: %.12f%s\n", i + 1, row_sums[i],
               (i == max_row) ? "  <- максимум" : "");
    }
    printf("q = норма alpha = максимальная сумма = %.12f\n\n", q);
    printf("Пояснение величины q:\n");
    printf("  q — норма матрицы alpha приведённой системы.\n");
    printf("  Условие q < 1 является достаточным условием сходимости и позволяет\n");
    printf("  использовать оценку погрешности <= q/(1-q) * max|dx|.\n");
    printf("  При q >= 1 эта оценка неприменима, но это не доказывает расходимость.\n");
    printf("  В этом случае сходимость контролируется по max|dx| и max|r_i|.\n\n");
}

static int gauss_seidel(const double a[N][N], const double b[N], double x[N],
                        double eps, int max_iter, int *iterations_out,
                        double *last_delta_out, double *residual_out,
                        double *error_estimate_out) {
    int iter = 0;
    int i = 0;
    int j = 0;
    double x_prev[N] = {0.0};
    double residual_vec[N] = {0.0};
    double q = calc_alpha_norm(a, NULL, NULL, NULL);

    for (i = 0; i < N; i++) {
        if (fabs(a[i][i]) < 1e-15) {
            return 0;
        }
        x[i] = 0.0;
    }

    for (iter = 1; iter <= max_iter; iter++) {
        double max_delta = 0.0;

        for (i = 0; i < N; i++) {
            x_prev[i] = x[i];
        }

        for (i = 0; i < N; i++) {
            double sum = b[i];
            for (j = 0; j < i; j++) {
                sum -= a[i][j] * x[j];
            }
            for (j = i + 1; j < N; j++) {
                sum -= a[i][j] * x_prev[j];
            }

            x[i] = sum / a[i][i];
            if (fabs(x[i] - x_prev[i]) > max_delta) {
                max_delta = fabs(x[i] - x_prev[i]);
            }
        }

        calc_residual(a, b, x, residual_vec, residual_out);
        *error_estimate_out = NAN;
        if (q < 1.0) {
            *error_estimate_out = q / (1.0 - q) * max_delta;
        }

        printf("Итерация %4d: [", iter);
        for (i = 0; i < N; i++) {
            printf("% .10f%s", x[i], (i + 1 < N) ? ", " : "");
        }
        printf("]\n");
        printf("  max|dx| = ");
        print_number_power10(max_delta, 12);
        printf(", max|r_i| = ");
        print_number_power10(*residual_out, 12);
        if (isfinite(*error_estimate_out)) {
            printf(", оценка погрешности <= ");
            print_number_power10(*error_estimate_out, 12);
            printf("\n");
        } else {
            printf(", оценка погрешности недоступна: q=%.6f >= 1\n", q);
        }

        if (((isfinite(*error_estimate_out) && *error_estimate_out < eps) ||
             (!isfinite(*error_estimate_out) && max_delta < eps)) &&
            *residual_out < eps) {
            *iterations_out = iter;
            *last_delta_out = max_delta;
            return 1;
        }
    }

    *iterations_out = max_iter;
    *last_delta_out = 0.0;
    calc_residual(a, b, x, residual_vec, residual_out);
    *error_estimate_out = NAN;
    return 0;
}

int run_task_6_2_1(void) {
    const double a[N][N] = {{2.67, 5.10, 3.31, 5.64, 4.76},
                            {4.44, 7.50, 4.67, 5.70, 6.14},
                            {5.33, 9.80, 8.67, 4.80, 7.33},
                            {3.56, 5.30, 4.15, 3.69, 3.25},
                            {1.78, 4.17, 2.67, 4.69, 3.75}};
    const double b[N] = {6.19, 6.95, 12.20, 5.97, 4.42};

    const int max_iter = 10000;
    int iterations = 0;
    int i = 0;
    double x[N] = {0.0};
    double residual_inf_norm = 0.0;
    double last_delta = 0.0;
    double error_estimate = NAN;
    double q_inf = 0.0;
    double alpha[N][N] = {{0.0}};
    double alpha_row_sums[N] = {0.0};
    int max_alpha_row = 0;
    double residual_vec[N] = {0.0};
    double eps = 1e-8;
    int converged = 0;

    print_system(a, b);

    printf("Метод Зейделя\n");
    printf("Введите eps (по умолчанию 1e-8): ");
    if (scanf("%lf", &eps) != 1 || eps <= 0.0) {
        eps = 1e-8;
        printf("Некорректный ввод, используется eps = ");
        print_number_power10(eps, 1);
        printf("\n");
    }

    if (!is_strict_diagonal_dominance(a)) {
        printf("Предупреждение: матрица не имеет строгого диагонального преобладания.\n");
        printf("Сходимость метода Зейделя этим условием не гарантируется.\n");
    }
    q_inf = calc_alpha_norm(a, alpha, alpha_row_sums, &max_alpha_row);
    print_alpha_details(a, b, alpha, alpha_row_sums, max_alpha_row, q_inf);
    if (q_inf < 1.0) {
        printf("Апостериорная оценка: q/(1-q)*||x(k)-x(k-1)||_inf\n");
    } else {
        printf("Апостериорная оценка по норме alpha неприменима, поскольку q >= 1.\n");
    }

    converged = gauss_seidel(a, b, x, eps, max_iter, &iterations, &last_delta,
                             &residual_inf_norm, &error_estimate);

    calc_residual(a, b, x, residual_vec, &residual_inf_norm);

    printf("\nРезультат:\n");
    if (converged) {
        printf("Метод сошёлся за %d итераций.\n", iterations);
    } else {
        printf("Метод не сошёлся за %d итераций.\n", iterations);
    }

    printf("Решение x:\n");
    for (i = 0; i < N; i++) {
        printf("x%d = %.12f\n", i + 1, x[i]);
    }

    printf("\nВектор невязки r = A*x - b:\n");
    for (i = 0; i < N; i++) {
        printf("r%d = ", i + 1);
        print_number_power10(residual_vec[i], 12);
        printf("\n");
    }

    printf("max|r_i| = ");
    print_number_power10(residual_inf_norm, 12);
    printf("\nПоследнее max|x(k)-x(k-1)| = ");
    print_number_power10(last_delta, 12);
    printf("\n");
    if (isfinite(error_estimate)) {
        printf("Теоретическая оценка погрешности <= ");
        print_number_power10(error_estimate, 12);
        printf("\n");
    } else {
        printf("Теоретическая оценка погрешности по q/(1-q) недоступна.\n");
    }

    return 0;
}
