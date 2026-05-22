#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "tasks/task_6_1/task_6_1.h"
#include "tasks/task_6_1/function_provider_6_1.h"
#include "tasks/task_6_1/solver_6_1.h"

#define MAX_EXPR_LEN 256
#define VARIANT_SPLIT_STEP 1.0

static double now_seconds(void) {
#ifdef _WIN32
    static LARGE_INTEGER freq;
    static int initialized = 0;
    LARGE_INTEGER counter;

    if (!initialized) {
        QueryPerformanceFrequency(&freq);
        initialized = 1;
    }
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)freq.QuadPart;
#else
    return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}

static int read_expression_line(const char *prompt, char *dst, size_t dst_size) {
    printf("%s", prompt);
    if (scanf(" %255[^\n]", dst) != 1) {
        return 0;
    }
    dst[dst_size - 1] = '\0';
    return 1;
}

int run_task_6_1(void) {
    double a = -5.0;
    double b = 5.0;
    const double split_step = VARIANT_SPLIT_STEP;
    double eps1 = 1e-8;
    double eps2 = 1e-8;
    int input_mode = 1;
    int function_mode = 1;
    char expr_f[MAX_EXPR_LEN];
    char expr_fp[MAX_EXPR_LEN];
    char expr_fpp[MAX_EXPR_LEN];
    RootInterval intervals[64];
    unsigned long split_f_calls = 0;
    int split_steps = 0;

#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif

    set_variant_function();

    printf("Комбинированный метод (Ньютона + хорд)\n");
    printf("Режим функции: 1 - вариант (sin(x)-0.2x), 2 - ручной ввод формул f, f', f'': ");
    if (scanf("%d", &function_mode) == 1 && function_mode == 2) {
        if (!read_expression_line("Введите f(x): ", expr_f, sizeof(expr_f)) ||
            !read_expression_line("Введите f'(x): ", expr_fp, sizeof(expr_fp)) ||
            !read_expression_line("Введите f''(x): ", expr_fpp, sizeof(expr_fpp)) ||
            !set_custom_functions(expr_f, expr_fp, expr_fpp)) {
            printf("Ошибка ввода формул. Используется функция варианта.\n");
            set_variant_function();
        }
    }

    printf("f(x): %s\n", get_function_description());
    printf("Режим ввода: 1 - стандартные параметры, 2 - ввод с консоли: ");
    if (scanf("%d", &input_mode) == 1 && input_mode == 2) {
        printf("Введите a b eps1 eps2: ");
        if (scanf("%lf %lf %lf %lf", &a, &b, &eps1, &eps2) != 4) {
            printf("Ошибка ввода. Используются стандартные параметры.\n");
            a = -5.0;
            b = 5.0;
            eps1 = 1e-8;
            eps2 = 1e-8;
        }
    }

    if (b <= a || eps1 <= 0.0 || eps2 <= 0.0) {
        printf("Некорректные параметры. Используются стандартные параметры.\n");
        a = -5.0;
        b = 5.0;
        eps1 = 1e-8;
        eps2 = 1e-8;
    }

    printf("Интервал: [%.3f, %.3f], шаг отделения: %.3f\n", a, b, split_step);
    printf("Точность: eps1 = %.1e, eps2 = %.1e\n\n", eps1, eps2);

    split_steps = (int)ceil((b - a) / split_step);
    printf("Количество шагов отделения корней: %d\n", split_steps);

    int roots_count = separate_roots(
        f, a, b, split_step, intervals,
        (int)(sizeof(intervals) / sizeof(intervals[0])), &split_f_calls);

    if (roots_count == 0) {
        printf("На интервале не найдено отрезков с корнями.\n");
        return 0;
    }

    printf("Найдено отрезков с корнями: %d\n", roots_count);
    printf("Вычисления f(x) на этапе отделения: %lu\n\n", split_f_calls);

    double total_time = 0.0;
    unsigned long total_f = split_f_calls;
    unsigned long total_fp = 0;
    unsigned long total_fpp = 0;

    for (int i = 0; i < roots_count; i++) {
        MethodStats stats;
        double root;
        double start = now_seconds();
        combined_method(f, f_prime, f_double_prime, intervals[i].left,
                        intervals[i].right, eps1, eps2, &root, &stats);
        double end = now_seconds();
        double elapsed = end - start;

        total_time += elapsed;
        total_f += stats.f_calls;
        total_fp += stats.fp_calls;
        total_fpp += stats.fpp_calls;

        printf("Корень #%d\n", i + 1);
        printf("  Отрезок: [%.10f, %.10f]\n", intervals[i].left, intervals[i].right);
        if (stats.converged) {
            printf("  ξ = %.12f\n", root);
            printf("  |f(ξ)| = %.3e\n", stats.residual_abs);
            printf("  Итерации n = %d\n", stats.iterations);
            printf("  Вычисления функций при уточнении: f=%lu, f'=%lu, f''=%lu\n",
                   stats.f_calls, stats.fp_calls, stats.fpp_calls);
            printf("  Время: %.9f сек\n", elapsed);
            if (isfinite(stats.last_dx)) {
                printf("  Последний шаг |x_n - x_(n-1)| = %.3e\n",
                       stats.last_dx);
            } else {
                printf("  Последний шаг |x_n - x_(n-1)| = не требуется\n");
            }
            printf("  Погрешность по аргументу <= %.3e\n", stats.argument_error);
            if (isfinite(stats.alpha)) {
                printf("  Параметр сходимости alpha = %.6f\n", stats.alpha);
            } else {
                printf("  Параметр сходимости alpha = н/д (мало данных)\n");
            }
        } else {
            printf("  Сходимость не достигнута.\n");
            printf("  Вызовы: f=%lu, f'=%lu, f''=%lu\n", stats.f_calls, stats.fp_calls,
                   stats.fpp_calls);
            printf("  Время: %.9f сек\n", elapsed);
        }
        printf("\n");
    }

    printf("ИТОГО\n");
    printf("  Шаги отделения корней: %d\n", split_steps);
    printf("  Суммарные вычисления функций: f=%lu, f'=%lu, f''=%lu\n", total_f, total_fp, total_fpp);
    printf("  Суммарное время: %.9f сек\n", total_time);
    return 0;
}
