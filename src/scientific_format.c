#include <math.h>
#include <stdio.h>
#include "scientific_format.h"

void print_number_power10(double value, int precision) {
    int exponent = 0;
    double mantissa = 0.0;

    if (value != 0.0) {
        exponent = (int)floor(log10(fabs(value)));
        mantissa = value / pow(10.0, exponent);
    }

    printf("% .*f · 10^(%d)", precision, mantissa, exponent);
}
