#include <math.h>
#include <stdio.h>
#include "tasks/task_6_2/task_6_2_1.h"

#define N 5

static void print_system(const double a[N][N], const double b[N]) {
    int i = 0;
    int j = 0;

    printf("System A*x=b (variant 9):\n");
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

static int gauss_seidel(const double a[N][N], const double b[N], double x[N],
                        double eps, int max_iter, int *iterations_out,
                        double *last_delta_out, double *residual_out) {
    int iter = 0;
    int i = 0;
    int j = 0;
    double x_prev[N] = {0.0};
    double residual_vec[N] = {0.0};

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

        printf("Iteration %4d: [", iter);
        for (i = 0; i < N; i++) {
            printf("% .10f%s", x[i], (i + 1 < N) ? ", " : "");
        }
        printf("]\n");

        calc_residual(a, b, x, residual_vec, residual_out);
        if (max_delta < eps && *residual_out < eps) {
            *iterations_out = iter;
            *last_delta_out = max_delta;
            return 1;
        }
    }

    *iterations_out = max_iter;
    *last_delta_out = 0.0;
    calc_residual(a, b, x, residual_vec, residual_out);
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
    double residual_vec[N] = {0.0};
    double eps = 1e-8;
    int converged = 0;

    print_system(a, b);

    printf("Seidel method\n");
    printf("Enter eps (default 1e-8): ");
    if (scanf("%lf", &eps) != 1 || eps <= 0.0) {
        eps = 1e-8;
        printf("Invalid input, using eps = %.1e\n", eps);
    }

    if (!is_strict_diagonal_dominance(a)) {
        printf("Warning: matrix has no strict diagonal dominance.\n");
        printf("Seidel convergence is not guaranteed.\n");
    }

    converged = gauss_seidel(a, b, x, eps, max_iter, &iterations, &last_delta,
                             &residual_inf_norm);

    calc_residual(a, b, x, residual_vec, &residual_inf_norm);

    printf("\nResult:\n");
    if (converged) {
        printf("Converged in %d iterations.\n", iterations);
    } else {
        printf("Did not converge in %d iterations.\n", iterations);
    }

    printf("Solution x:\n");
    for (i = 0; i < N; i++) {
        printf("x%d = %.12f\n", i + 1, x[i]);
    }

    printf("\nResidual vector r = A*x - b:\n");
    for (i = 0; i < N; i++) {
        printf("r%d = %.12e\n", i + 1, residual_vec[i]);
    }

    printf("||r||_inf = %.12e\n", residual_inf_norm);
    printf("Last max|x(k)-x(k-1)| = %.12e\n", last_delta);

    return 0;
}
