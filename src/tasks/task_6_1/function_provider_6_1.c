#include <math.h>
#include <string.h>
#include "tasks/task_6_1/expression_6_1.h"
#include "tasks/task_6_1/function_provider_6_1.h"

#define MAX_EXPR_LEN 256

typedef enum {
    MODE_VARIANT = 1,
    MODE_EXPR = 2
} FunctionMode;

static FunctionMode g_mode = MODE_VARIANT;
static char g_expr_f[MAX_EXPR_LEN] = "sin(x)-0.2*x";
static char g_expr_fp[MAX_EXPR_LEN] = "cos(x)-0.2";
static char g_expr_fpp[MAX_EXPR_LEN] = "-sin(x)";
static char g_description[MAX_EXPR_LEN] = "sin(x) - 0.2x";

void set_variant_function(void) {
    g_mode = MODE_VARIANT;
    strcpy(g_expr_f, "sin(x)-0.2*x");
    strcpy(g_expr_fp, "cos(x)-0.2");
    strcpy(g_expr_fpp, "-sin(x)");
    strcpy(g_description, "sin(x) - 0.2x");
}

int set_custom_functions(const char *fx, const char *fpx, const char *fppx) {
    double test_value = 0.0;

    if (fx == NULL || fpx == NULL || fppx == NULL) {
        return 0;
    }
    if (!eval_expression(fx, 0.0, &test_value)) {
        return 0;
    }
    if (!eval_expression(fpx, 0.0, &test_value)) {
        return 0;
    }
    if (!eval_expression(fppx, 0.0, &test_value)) {
        return 0;
    }

    strncpy(g_expr_f, fx, MAX_EXPR_LEN - 1);
    g_expr_f[MAX_EXPR_LEN - 1] = '\0';
    strncpy(g_expr_fp, fpx, MAX_EXPR_LEN - 1);
    g_expr_fp[MAX_EXPR_LEN - 1] = '\0';
    strncpy(g_expr_fpp, fppx, MAX_EXPR_LEN - 1);
    g_expr_fpp[MAX_EXPR_LEN - 1] = '\0';

    strncpy(g_description, g_expr_f, MAX_EXPR_LEN - 1);
    g_description[MAX_EXPR_LEN - 1] = '\0';
    g_mode = MODE_EXPR;
    return 1;
}

const char *get_function_description(void) {
    return g_description;
}

double f(double x) {
    if (g_mode == MODE_EXPR) {
        double value = NAN;
        if (!eval_expression(g_expr_f, x, &value)) {
            return NAN;
        }
        return value;
    }
    return sin(x) - 0.2 * x;
}

double f_prime(double x) {
    if (g_mode == MODE_EXPR) {
        double value = NAN;
        if (!eval_expression(g_expr_fp, x, &value)) {
            return NAN;
        }
        return value;
    }
    return cos(x) - 0.2;
}

double f_double_prime(double x) {
    if (g_mode == MODE_EXPR) {
        double value = NAN;
        if (!eval_expression(g_expr_fpp, x, &value)) {
            return NAN;
        }
        return value;
    }
    return -sin(x);
}
