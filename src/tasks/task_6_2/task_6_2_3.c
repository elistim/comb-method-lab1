#include <math.h>
#include <stdio.h>
#include "scientific_format.h"
#include "tasks/task_6_2/khaletsky.h"
#include "tasks/task_6_2/task_6_2_3.h"

static void invert_from_khaletsky(
    const double l[KHALETSKY_N][KHALETSKY_N],
    const double u[KHALETSKY_N][KHALETSKY_N],
    double inv_a[KHALETSKY_N][KHALETSKY_N]) {
    int col = 0;
    int i = 0;
    int k = 0;
    double y[KHALETSKY_N];
    double x[KHALETSKY_N];

    for (i = 0; i < KHALETSKY_N; i++) {
        for (k = 0; k < KHALETSKY_N; k++) {
            inv_a[i][k] = 0.0;
        }
    }

    for (col = 0; col < KHALETSKY_N; col++) {
        for (i = 0; i < KHALETSKY_N; i++) {
            double sum = (i == col) ? 1.0 : 0.0;
            for (k = 0; k < i; k++) {
                sum -= l[i][k] * y[k];
            }
            y[i] = sum / l[i][i];
        }

        for (i = KHALETSKY_N - 1; i >= 0; i--) {
            double sum = y[i];
            for (k = i + 1; k < KHALETSKY_N; k++) {
                sum -= u[i][k] * x[k];
            }
            x[i] = sum;
        }

        for (i = 0; i < KHALETSKY_N; i++) {
            inv_a[i][col] = x[i];
        }

        {
            char title_y[64];
            char title_x[64];
            snprintf(title_y, sizeof(title_y), "Столбец %d: y из L*y=e:",
                     col + 1);
            snprintf(title_x, sizeof(title_x), "Столбец %d: x из U*x=y:",
                     col + 1);
            print_vector(title_y, y);
            print_vector(title_x, x);
        }
    }
    printf("\n");
}

static double residual_norm_inf(
    const double a[KHALETSKY_N][KHALETSKY_N],
    const double inv_a[KHALETSKY_N][KHALETSKY_N],
    double r[KHALETSKY_N][KHALETSKY_N]) {
    double aa_inv[KHALETSKY_N][KHALETSKY_N];
    double max_row_sum = 0.0;
    int i = 0;
    int j = 0;

    multiply_matrix(a, inv_a, aa_inv);

    for (i = 0; i < KHALETSKY_N; i++) {
        double row_sum = 0.0;
        for (j = 0; j < KHALETSKY_N; j++) {
            double ij = (i == j) ? 1.0 : 0.0;
            r[i][j] = aa_inv[i][j] - ij;
            row_sum += fabs(r[i][j]);
        }
        if (row_sum > max_row_sum) {
            max_row_sum = row_sum;
        }
    }

    return max_row_sum;
}

int run_task_6_2_3(void) {
    const double a[KHALETSKY_N][KHALETSKY_N] = {
        {2.67, 5.10, 3.31, 5.64, 4.76},
        {4.44, 7.50, 4.67, 5.70, 6.14},
        {5.33, 9.80, 8.67, 4.80, 7.33},
        {3.56, 5.30, 4.15, 3.69, 3.25},
        {1.78, 4.17, 2.67, 4.69, 3.75}};
    double l[KHALETSKY_N][KHALETSKY_N];
    double u[KHALETSKY_N][KHALETSKY_N];
    double inv_a[KHALETSKY_N][KHALETSKY_N];
    double r[KHALETSKY_N][KHALETSKY_N];
    double norm_inf = 0.0;

    printf("6.2.3 Вычисление обратной матрицы методом Халецкого\n");
    print_matrix("A =", a);

    if (!khaletsky_decompose(a, l, u, print_khaletsky_step)) {
        printf("Не удалось выполнить разложение Халецкого.\n");
        return 1;
    }

    print_matrix("Итоговая L =", l);
    print_matrix("Итоговая U =", u);

    invert_from_khaletsky(l, u, inv_a);
    print_matrix("A^-1 =", inv_a);

    norm_inf = residual_norm_inf(a, inv_a, r);
    print_matrix_scientific("Погрешность R = A * A^-1 - I =", r);
    printf("max row sum = ");
    print_number_power10(norm_inf, 12);
    printf("\n");

    return 0;
}
