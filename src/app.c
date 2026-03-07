#include <stdio.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#endif
#include "app.h"
#include "tasks/task_6_1/task_6_1.h"
#include "tasks/task_6_2/task_6_2_1.h"
#include "tasks/task_6_2/task_6_2_2.h"
#include "tasks/task_6_2/task_6_2_3.h"

static int run_selected_task(const char *task_id) {
    if (strcmp(task_id, "6.1") == 0) {
        return run_task_6_1();
    }
    if (strcmp(task_id, "6.2.1") == 0) {
        return run_task_6_2_1();
    }
    if (strcmp(task_id, "6.2.2") == 0) {
        return run_task_6_2_2();
    }
    if (strcmp(task_id, "6.2.3") == 0) {
        return run_task_6_2_3();
    }

    printf("Неизвестный номер задания: %s\n", task_id);
    return 1;
}

int run_application(void) {
    char task_id[16] = {0};

#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
#endif

    printf("Выберите задание:\n");
    printf("  6.1   - Решение уравнений с одной переменной\n");
    printf("  6.2.1 - Решение систем линейных уравнений\n");
    printf("  6.2.2 - Вычисление определителей матриц\n");
    printf("  6.2.3 - Вычисление обратной матрицы\n");
    printf("Введите номер задания: ");

    if (scanf("%15s", task_id) != 1) {
        printf("Ошибка чтения номера задания.\n");
        return 1;
    }

    return run_selected_task(task_id);
}
