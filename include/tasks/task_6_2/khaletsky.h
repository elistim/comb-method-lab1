#ifndef KHALETSKY_H
#define KHALETSKY_H

#define KHALETSKY_N 5

typedef void (*KhaletskyStepPrinter)(int step,
                                     const double l[KHALETSKY_N][KHALETSKY_N],
                                     const double u[KHALETSKY_N][KHALETSKY_N]);

void print_matrix(const char *title,
                  const double m[KHALETSKY_N][KHALETSKY_N]);
void print_matrix_scientific(const char *title,
                             const double m[KHALETSKY_N][KHALETSKY_N]);
void print_vector(const char *title, const double v[KHALETSKY_N]);
void multiply_matrix(const double a[KHALETSKY_N][KHALETSKY_N],
                     const double b[KHALETSKY_N][KHALETSKY_N],
                     double c[KHALETSKY_N][KHALETSKY_N]);
int khaletsky_decompose(const double a[KHALETSKY_N][KHALETSKY_N],
                        double l[KHALETSKY_N][KHALETSKY_N],
                        double u[KHALETSKY_N][KHALETSKY_N],
                        KhaletskyStepPrinter step_printer);
double khaletsky_det(const double l[KHALETSKY_N][KHALETSKY_N]);
void print_khaletsky_step(int step,
                          const double l[KHALETSKY_N][KHALETSKY_N],
                          const double u[KHALETSKY_N][KHALETSKY_N]);

#endif // KHALETSKY_H
