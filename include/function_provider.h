#ifndef FUNCTION_PROVIDER_H
#define FUNCTION_PROVIDER_H

void set_variant_function(void);
int set_custom_functions(const char *fx, const char *fpx, const char *fppx);
const char *get_function_description(void);

double f(double x);
double f_prime(double x);
double f_double_prime(double x);

#endif // FUNCTION_PROVIDER_H
