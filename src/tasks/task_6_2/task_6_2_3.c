#include <math.h>
#include <stdio.h>
#include "tasks/task_6_2/task_6_2_3.h"

#define N 5
#define EPS 1e-12

static void print_matrix(const char *title, const double m[N][N]) {
    int i = 0;
    int j = 0;

    printf("%s\n", title);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%13.8f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

static void transpose(const double a[N][N], double at[N][N]) {
    int i = 0;
    int j = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            at[j][i] = a[i][j];
        }
    }
}

static void multiply(const double a[N][N], const double b[N][N], double c[N][N]) {
    int i = 0;
    int j = 0;
    int k = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            double s = 0.0;
            for (k = 0; k < N; k++) {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
}

static int cholesky_decompose(const double c[N][N], double l[N][N]) {
    int i = 0;
    int j = 0;
    int k = 0;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            l[i][j] = 0.0;
        }
    }

    for (i = 0; i < N; i++) {
        for (j = 0; j <= i; j++) {
            double sum = c[i][j];

            for (k = 0; k < j; k++) {
                sum -= l[i][k] * l[j][k];
            }

            if (i == j) {
                if (sum <= EPS) {
                    return 0;
                }
                l[i][j] = sqrt(sum);
            } else {
                if (fabs(l[j][j]) < EPS) {
                    return 0;
                }
                l[i][j] = sum / l[j][j];
            }
        }
    }

    return 1;
}

static void invert_spd_from_cholesky(const double l[N][N], double inv_c[N][N]) {
    int col = 0;
    int i = 0;
    int k = 0;
    double y[N];
    double x[N];

    for (i = 0; i < N; i++) {
        for (k = 0; k < N; k++) {
            inv_c[i][k] = 0.0;
        }
    }

    for (col = 0; col < N; col++) {
        for (i = 0; i < N; i++) {
            double sum = (i == col) ? 1.0 : 0.0;
            for (k = 0; k < i; k++) {
                sum -= l[i][k] * y[k];
            }
            y[i] = sum / l[i][i];
        }

        for (i = N - 1; i >= 0; i--) {
            double sum = y[i];
            for (k = i + 1; k < N; k++) {
                sum -= l[k][i] * x[k];
            }
            x[i] = sum / l[i][i];
        }

        for (i = 0; i < N; i++) {
            inv_c[i][col] = x[i];
        }
    }
}

static double residual_norm_inf(const double a[N][N], const double inv_a[N][N],
                                double r[N][N]) {
    double aa_inv[N][N];
    double max_row_sum = 0.0;
    int i = 0;
    int j = 0;

    multiply(a, inv_a, aa_inv);

    for (i = 0; i < N; i++) {
        double row_sum = 0.0;
        for (j = 0; j < N; j++) {
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
    const double a[N][N] = {{2.67, 5.10, 3.31, 5.64, 4.76},
                            {4.44, 7.50, 4.67, 5.70, 6.14},
                            {5.33, 9.80, 8.67, 4.80, 7.33},
                            {3.56, 5.30, 4.15, 3.69, 3.25},
                            {1.78, 4.17, 2.67, 4.69, 3.75}};
    double at[N][N];
    double c[N][N];
    double l[N][N];
    double inv_c[N][N];
    double inv_a[N][N];
    double r[N][N];
    double norm_inf = 0.0;

    printf("6.2.3 Вычисление обратной матрицы методом Халецкого\n");
    print_matrix("A =", a);

    transpose(a, at);
    print_matrix("A^T =", at);

    multiply(at, a, c);
    print_matrix("C = A^T * A =", c);

    if (!cholesky_decompose(c, l)) {
        printf("Не удалось выполнить разложение Халецкого (C не SPD).\n");
        return 1;
    }

    print_matrix("L (C = L * L^T) =", l);

    invert_spd_from_cholesky(l, inv_c);
    print_matrix("(A^T * A)^-1 =", inv_c);

    multiply(inv_c, at, inv_a);
    print_matrix("A^-1 = (A^T * A)^-1 * A^T =", inv_a);

    norm_inf = residual_norm_inf(a, inv_a, r);
    print_matrix("Невязка R = A * A^-1 - I =", r);
    printf("||R||_inf = %.12e\n", norm_inf);

    return 0;
}
