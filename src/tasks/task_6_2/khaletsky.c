#include <math.h>
#include <stdio.h>
#include "tasks/task_6_2/khaletsky.h"

#define EPS_PIVOT 1e-12

void print_matrix(const char *title,
                  const double m[KHALETSKY_N][KHALETSKY_N]) {
    int i = 0;
    int j = 0;

    printf("%s\n", title);
    for (i = 0; i < KHALETSKY_N; i++) {
        for (j = 0; j < KHALETSKY_N; j++) {
            printf("%13.8f ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_scientific(const char *title,
                             const double m[KHALETSKY_N][KHALETSKY_N]) {
    int i = 0;
    int j = 0;

    printf("%s\n", title);
    for (i = 0; i < KHALETSKY_N; i++) {
        for (j = 0; j < KHALETSKY_N; j++) {
            printf("%13.6e ", m[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vector(const char *title, const double v[KHALETSKY_N]) {
    int i = 0;

    printf("%s [", title);
    for (i = 0; i < KHALETSKY_N; i++) {
        printf("% .12f%s", v[i], (i + 1 < KHALETSKY_N) ? ", " : "");
    }
    printf("]\n");
}

void multiply_matrix(const double a[KHALETSKY_N][KHALETSKY_N],
                     const double b[KHALETSKY_N][KHALETSKY_N],
                     double c[KHALETSKY_N][KHALETSKY_N]) {
    int i = 0;
    int j = 0;
    int k = 0;

    for (i = 0; i < KHALETSKY_N; i++) {
        for (j = 0; j < KHALETSKY_N; j++) {
            double s = 0.0;
            for (k = 0; k < KHALETSKY_N; k++) {
                s += a[i][k] * b[k][j];
            }
            c[i][j] = s;
        }
    }
}

int khaletsky_decompose(const double a[KHALETSKY_N][KHALETSKY_N],
                        double l[KHALETSKY_N][KHALETSKY_N],
                        double u[KHALETSKY_N][KHALETSKY_N],
                        KhaletskyStepPrinter step_printer) {
    int i = 0;
    int j = 0;
    int k = 0;

    for (i = 0; i < KHALETSKY_N; i++) {
        for (j = 0; j < KHALETSKY_N; j++) {
            l[i][j] = 0.0;
            u[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (k = 0; k < KHALETSKY_N; k++) {
        for (i = k; i < KHALETSKY_N; i++) {
            double sum = 0.0;
            for (j = 0; j < k; j++) {
                sum += l[i][j] * u[j][k];
            }
            l[i][k] = a[i][k] - sum;
        }

        if (fabs(l[k][k]) < EPS_PIVOT) {
            return 0;
        }

        for (j = k + 1; j < KHALETSKY_N; j++) {
            double sum = 0.0;
            for (i = 0; i < k; i++) {
                sum += l[k][i] * u[i][j];
            }
            u[k][j] = (a[k][j] - sum) / l[k][k];
        }

        if (step_printer != NULL) {
            step_printer(k + 1, l, u);
        }
    }

    return 1;
}

double khaletsky_det(const double l[KHALETSKY_N][KHALETSKY_N]) {
    double det = 1.0;
    int i = 0;

    for (i = 0; i < KHALETSKY_N; i++) {
        det *= l[i][i];
    }

    return det;
}

void print_khaletsky_step(int step,
                          const double l[KHALETSKY_N][KHALETSKY_N],
                          const double u[KHALETSKY_N][KHALETSKY_N]) {
    char title_l[64];
    char title_u[64];

    snprintf(title_l, sizeof(title_l), "L после шага k=%d:", step);
    snprintf(title_u, sizeof(title_u), "U после шага k=%d:", step);
    print_matrix(title_l, l);
    print_matrix(title_u, u);
}
