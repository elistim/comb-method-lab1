#include <math.h>
#include <stdio.h>
#include "tasks/task_6_2/task_6_2_2.h"

#define N 5
#define EPS_PIVOT 1e-12

static void copy_matrix(const double src[N][N], double dst[N][N]) {
    int i = 0;
    int j = 0;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            dst[i][j] = src[i][j];
        }
    }
}

static void print_matrix(const char *title, const double m[N][N]) {
    int i = 0;
    int j = 0;
    printf("%s\n", title);
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%12.6f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

static void swap_rows(double m[N][N], int r1, int r2) {
    int j = 0;
    for (j = 0; j < N; j++) {
        double t = m[r1][j];
        m[r1][j] = m[r2][j];
        m[r2][j] = t;
    }
}

static int det_by_gauss(const double a[N][N], double *det_out) {
    double u[N][N];
    int sign = 1;
    int i = 0;
    int j = 0;
    int k = 0;

    copy_matrix(a, u);
    print_matrix("Исходная матрица A:", u);

    for (k = 0; k < N; k++) {
        int pivot_row = k;
        double max_abs = fabs(u[k][k]);

        for (i = k + 1; i < N; i++) {
            double v = fabs(u[i][k]);
            if (v > max_abs) {
                max_abs = v;
                pivot_row = i;
            }
        }

        if (max_abs < EPS_PIVOT) {
            *det_out = 0.0;
            return 0;
        }

        if (pivot_row != k) {
            swap_rows(u, pivot_row, k);
            sign = -sign;
            printf("Перестановка строк %d <-> %d\n", pivot_row + 1, k + 1);
        }

        for (i = k + 1; i < N; i++) {
            double factor = u[i][k] / u[k][k];
            for (j = k; j < N; j++) {
                u[i][j] -= factor * u[k][j];
            }
        }

        {
            char title[64];
            snprintf(title, sizeof(title), "После шага Гаусса k=%d:", k + 1);
            print_matrix(title, u);
        }
    }

    *det_out = (double)sign;
    for (i = 0; i < N; i++) {
        *det_out *= u[i][i];
    }
    return 1;
}

int run_task_6_2_2(void) {
    const double a[N][N] = {{2.67, 5.10, 3.31, 5.64, 4.76},
                            {4.44, 7.50, 4.67, 5.70, 6.14},
                            {5.33, 9.80, 8.67, 4.80, 7.33},
                            {3.56, 5.30, 4.15, 3.69, 3.25},
                            {1.78, 4.17, 2.67, 4.69, 3.75}};
    double det_gauss = 0.0;

    printf("6.2.2 Вычисление определителя матрицы методом Гаусса\n");
    printf("Матрица взята из задания 6.2.1 (вариант 9)\n\n");

    if (det_by_gauss(a, &det_gauss)) {
        printf("det(A) методом Гаусса = %.12f\n", det_gauss);
    } else {
        printf("Матрица вырождена, det(A)=0\n");
    }

    return 0;
}
