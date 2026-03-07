#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tasks/task_6_1/expression_6_1.h"

typedef struct {
    const char *s;
    double x;
    int error;
} ExprParser;

static void skip_spaces(ExprParser *p) {
    while (*p->s != '\0' && isspace((unsigned char)*p->s)) {
        p->s++;
    }
}

static int consume_word(ExprParser *p, const char *word) {
    size_t len = strlen(word);
    if (strncmp(p->s, word, len) == 0) {
        p->s += len;
        return 1;
    }
    return 0;
}

static double parse_expression(ExprParser *p);

static double parse_primary(ExprParser *p) {
    skip_spaces(p);

    if (*p->s == '(') {
        p->s++;
        double v = parse_expression(p);
        skip_spaces(p);
        if (*p->s != ')') {
            p->error = 1;
            return 0.0;
        }
        p->s++;
        return v;
    }

    if (*p->s == 'x' || *p->s == 'X') {
        p->s++;
        return p->x;
    }

    if (isalpha((unsigned char)*p->s)) {
        double arg = 0.0;

        if (consume_word(p, "sin")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return sin(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "cos")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return cos(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "tan")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return tan(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "exp")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return exp(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "log")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return log(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "sqrt")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return sqrt(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }
        if (consume_word(p, "abs")) {
            skip_spaces(p);
            if (*p->s == '(') {
                p->s++;
                arg = parse_expression(p);
                skip_spaces(p);
                if (*p->s == ')') {
                    p->s++;
                    return fabs(arg);
                }
            }
            p->error = 1;
            return 0.0;
        }

        p->error = 1;
        return 0.0;
    }

    char *end_ptr = NULL;
    double v = strtod(p->s, &end_ptr);
    if (end_ptr == p->s) {
        p->error = 1;
        return 0.0;
    }
    p->s = end_ptr;
    return v;
}

static double parse_unary(ExprParser *p) {
    skip_spaces(p);
    if (*p->s == '+') {
        p->s++;
        return parse_unary(p);
    }
    if (*p->s == '-') {
        p->s++;
        return -parse_unary(p);
    }
    return parse_primary(p);
}

static double parse_power(ExprParser *p) {
    double left = parse_unary(p);
    skip_spaces(p);
    if (*p->s == '^') {
        p->s++;
        return pow(left, parse_power(p));
    }
    return left;
}

static double parse_term(ExprParser *p) {
    double v = parse_power(p);
    while (!p->error) {
        skip_spaces(p);
        if (*p->s == '*') {
            p->s++;
            v *= parse_power(p);
        } else if (*p->s == '/') {
            p->s++;
            v /= parse_power(p);
        } else {
            break;
        }
    }
    return v;
}

static double parse_expression(ExprParser *p) {
    double v = parse_term(p);
    while (!p->error) {
        skip_spaces(p);
        if (*p->s == '+') {
            p->s++;
            v += parse_term(p);
        } else if (*p->s == '-') {
            p->s++;
            v -= parse_term(p);
        } else {
            break;
        }
    }
    return v;
}

int eval_expression(const char *expr, double x, double *value_out) {
    ExprParser p;
    if (expr == NULL || value_out == NULL) {
        return 0;
    }

    p.s = expr;
    p.x = x;
    p.error = 0;

    double v = parse_expression(&p);
    skip_spaces(&p);

    if (p.error || *p.s != '\0' || !isfinite(v)) {
        return 0;
    }

    *value_out = v;
    return 1;
}
