#ifndef SOLVER_6_1_H
#define SOLVER_6_1_H

typedef struct {
    double left;
    double right;
} RootInterval;

typedef struct {
    int iterations;
    unsigned long f_calls;
    unsigned long fp_calls;
    unsigned long fpp_calls;
    double alpha;
    double last_dx;
    int converged;
} MethodStats;

int separate_roots(double (*func)(double), double a, double b, double step,
                   RootInterval *intervals, int max_count,
                   unsigned long *f_calls);

/* Комбинированный метод (Ньютона + хорд) для задания 6.1. */
void combined_method(double (*func)(double), double (*func_prime)(double),
                     double (*func_double_prime)(double), double a, double b,
                     double eps1, double eps2, double *root,
                     MethodStats *stats);

#endif // SOLVER_6_1_H
