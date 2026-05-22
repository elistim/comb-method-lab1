#include <stdio.h>
#include "tasks/task_6_2/khaletsky.h"
#include "tasks/task_6_2/task_6_2_2.h"

static int det_by_khaletsky(const double a[KHALETSKY_N][KHALETSKY_N],
                            double *det_out) {
    double l[KHALETSKY_N][KHALETSKY_N];
    double u[KHALETSKY_N][KHALETSKY_N];

    print_matrix("Исходная матрица A:", a);
    if (!khaletsky_decompose(a, l, u, print_khaletsky_step)) {
        *det_out = 0.0;
        return 0;
    }

    *det_out = khaletsky_det(l);
    return 1;
}

int run_task_6_2_2(void) {
    const double a[KHALETSKY_N][KHALETSKY_N] = {
        {2.67, 5.10, 3.31, 5.64, 4.76},
        {4.44, 7.50, 4.67, 5.70, 6.14},
        {5.33, 9.80, 8.67, 4.80, 7.33},
        {3.56, 5.30, 4.15, 3.69, 3.25},
        {1.78, 4.17, 2.67, 4.69, 3.75}};
    double det = 0.0;

    printf("6.2.2 Вычисление определителя матрицы методом Халецкого\n");
    printf("Матрица взята из задания 6.2.1 (вариант 9)\n\n");

    if (det_by_khaletsky(a, &det)) {
        printf("det(A) методом Халецкого = %.12f\n", det);
    } else {
        printf("Матрица вырождена, det(A)=0\n");
    }

    return 0;
}
