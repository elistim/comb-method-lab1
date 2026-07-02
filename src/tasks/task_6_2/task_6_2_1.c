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

static double calc_seidel_q_inf(const double a[N][N], double b_out[N][N],
                                double row_sums_out[N], int *max_row_out) {
    double b_iter[N][N] = {{0.0}};
    double q = 0.0;
    int max_row = 0;
    int col = 0;
    int i = 0;
    int j = 0;

    for (col = 0; col < N; col++) {
        for (i = 0; i < N; i++) {
            double sum = (col > i) ? -a[i][col] : 0.0;

            for (j = 0; j < i; j++) {
                sum -= a[i][j] * b_iter[j][col];
            }

            b_iter[i][col] = sum / a[i][i];
        }
    }

    for (i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (j = 0; j < N; j++) {
            row_sum += fabs(b_iter[i][j]);
            if (b_out != NULL) {
                b_out[i][j] = b_iter[i][j];
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

static void print_seidel_matrix_details(const double a[N][N],
                                        const double b_iter[N][N],
                                        const double row_sums[N],
                                        int max_row, double q) {
    int i = 0;
    int j = 0;

    printf("\nРазложение A = (D + L) + U, где:\n");
    printf("D + L — диагональная и нижняя треугольная части:\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < N; j++) {
            printf("% .6f%s", (j <= i) ? a[i][j] : 0.0,
                   (j + 1 < N) ? "  " : "");
        }
        printf("]\n");
    }

    printf("U — строго верхняя треугольная часть:\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < N; j++) {
            printf("% .6f%s", (j > i) ? a[i][j] : 0.0,
                   (j + 1 < N) ? "  " : "");
        }
        printf("]\n");
    }

    printf("Для каждого столбца j решается система (D + L) * B_j = -U_j.\n");
    printf("Системы для вычисления столбцов B_Seidel:\n");
    for (j = 0; j < N; j++) {
        int row = 0;
        int k = 0;
        printf("\nСтолбец B_%d:\n", j + 1);
        for (row = 0; row < N; row++) {
            printf("  ");
            for (k = 0; k <= row; k++) {
                if (k > 0) {
                    printf(" + ");
                }
                printf("%.2f*B_%d%d", a[row][k], k + 1, j + 1);
            }
            printf(" = %.2f\n", (j > row) ? -a[row][j] : 0.0);
        }
        printf("  Решение столбца: ");
        for (row = 0; row < N; row++) {
            printf("B_%d%d=%.6f%s", row + 1, j + 1, b_iter[row][j],
                   (row + 1 < N) ? ", " : "\n");
        }
    }

    printf("\nПолученные столбцы образуют B_Seidel = -(D + L)^(-1) * U:\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < N; j++) {
            printf("% .12f%s", b_iter[i][j], (j + 1 < N) ? "  " : "");
        }
        printf("]\n");
    }

    printf("\nСуммы модулей элементов строк B_Seidel:\n");
    for (i = 0; i < N; i++) {
        printf("s_%d = ", i + 1);
        for (j = 0; j < N; j++) {
            printf("|%.12f|%s", b_iter[i][j], (j + 1 < N) ? " + " : "");
        }
        printf(" = %.12f%s\n", row_sums[i],
               (i == max_row) ? "  <- максимум" : "");
    }
    printf("q = ||B_Seidel||_inf = максимум(s_1, ..., s_%d) = %.12f\n\n", N, q);
    printf("Пояснение величины q:\n");
    printf("  q — бесконечная норма итерационной матрицы метода Зейделя.\n");
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
    double q = calc_seidel_q_inf(a, NULL, NULL, NULL);

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
    double b_seidel[N][N] = {{0.0}};
    double b_row_sums[N] = {0.0};
    int max_b_row = 0;
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
    q_inf = calc_seidel_q_inf(a, b_seidel, b_row_sums, &max_b_row);
    print_seidel_matrix_details(a, b_seidel, b_row_sums, max_b_row, q_inf);
    if (q_inf < 1.0) {
        printf("Апостериорная оценка: q/(1-q)*||x(k)-x(k-1)||_inf\n");
    } else {
        printf("Апостериорная оценка по ||B||_inf неприменима, поскольку q >= 1.\n");
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
