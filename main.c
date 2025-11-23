#define _GNU_SOURCE
#include "main.h"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/resource.h>

static StateMath STATE_MATH = START;

complex_t cmul(complex_t a, complex_t b)
{
    complex_t r;
    r.re = a.re * b.re - a.im * b.im;
    r.im = a.re * b.im + a.im * b.re;
    return r;
}

// комплексное сложение
complex_t cadd(complex_t a, complex_t b)
{
    complex_t r;
    r.re = a.re + b.re;
    r.im = a.im + b.im;
    return r;
}

// Комплексное сопряжение
complex_t cconj(complex_t a)
{
    complex_t r;
    r.re = a.re;
    r.im = -a.im;
    return r;
}

// Переменные
double dx = 0, dy = 0;
int Nx = 0, Ny = 0;
double Tx = 0, Ty = 0;
double Tx_min = 0, Tx_max = 0, Ty_min = 0, Ty_max = 0;
double T_start = 0, T_step = 0, T_end = 0;
double TimeSimul[255];
char SNR[32] = "NONE";
int N_Signal = 0;
double F[64]; // массивы с запасом (до 64 сигналов)
double x[64];
double y[64];
double A[64];

#define TIME_NUMBER_STEPS (int)((T_end - T_start) / T_step)

int T = 0; // число временных отсчётов
int M = 0;

int lenPhi = 0;
int lenTheta = 0;


double read_cpu_temp() {
    FILE *fp = fopen("/sys/class/thermal/thermal_zone0/temp", "r");
    if (!fp) return -1.0;

    int temp_milli;
    if (fscanf(fp, "%d", &temp_milli) != 1) {
        fclose(fp);
        return -1.0;
    }
    fclose(fp);
    return temp_milli / 1000.0;  // в °C
}

/**
 * @brief Точка входа в программу.
 *
 * @details
 * Функция инициализирует основные параметры, запускает вычисления
 * и управляет основным циклом программы.
 * В зависимости от переданных аргументов командной строки
 * может изменять поведение (например, задавать размеры сетки, количество потоков и т.п.).
 *
 * @param argc Количество аргументов командной строки.
 * @param argv Массив строковых аргументов командной строки.
 *
 * @return Код завершения программы:
 * - 0 — успешное завершение;
 * - ненулевое значение — возникла ошибка.
 *
 * @note Это стандартная функция `main()`, которая вызывается при старте программы.
 */
int main(int argc, char *argv[])
{

#ifndef OFF_ARGS

    if (argc != 2)
    {
        fprintf(stderr, "Ошибка: ожидается ровно один аргумент (1D или 2D).\n");
        fprintf(stderr, "Пример: %s 1D\n", argv[0]);
        return EXIT_FAILURE;
    }
#endif
    FILE *fp = fopen("Beam.params", "r");
    if (!fp)
    {
        perror("Ошибка открытия Beam.params");
        return EXIT_FAILURE;
    }

struct sched_param param;
    param.sched_priority = 90; // приоритет 1–99 для RT

    if (sched_setscheduler(0, SCHED_FIFO, &param) == -1) {
        perror("sched_setscheduler");
    } else {
        printf("Real-time priority set: %d\n", param.sched_priority);
    }

    char line[MAX_LINE];
    while (fgets(line, sizeof(line), fp))
    {
        char key[64], value[64];

        if (sscanf(line, "%63s = %63s", key, value) == 2)
        {
            if (strcmp(key, "dx") == 0)
                dx = atof(value);
            else if (strcmp(key, "dy") == 0)
                dy = atof(value);
            else if (strcmp(key, "Nx") == 0)
                Nx = atoi(value);
            else if (strcmp(key, "Ny") == 0)
                Ny = atoi(value);

            else if (strcmp(key, "Tx") == 0)
                Tx = atof(value);
            else if (strcmp(key, "Ty") == 0)
                Ty = atof(value);
            else if (strcmp(key, "Tx_min") == 0)
                Tx_min = atof(value);
            else if (strcmp(key, "Tx_max") == 0)
                Tx_max = atof(value);
            else if (strcmp(key, "Ty_min") == 0)
                Ty_min = atof(value);
            else if (strcmp(key, "Ty_max") == 0)
                Ty_max = atof(value);

            else if (strcmp(key, "T_start") == 0)
                T_start = atof(value);
            else if (strcmp(key, "T_step") == 0)
                T_step = atof(value);
            else if (strcmp(key, "T_end") == 0)
                T_end = atof(value);

            else if (strcmp(key, "SNR") == 0)
                strncpy(SNR, value, sizeof(SNR));

            else if (strcmp(key, "N_Signal") == 0)
            {
                N_Signal = atoi(value);
            }
            else
            {
                // читаем параметры сигналов в цикле
                for (int i = 0; i < N_Signal; i++)
                {
                    char keyF[16], keyX[16], keyY[16], keyA[16];
                    sprintf(keyF, "F%d", i + 1);
                    sprintf(keyX, "x%d", i + 1);
                    sprintf(keyY, "y%d", i + 1);
                    sprintf(keyA, "A%d", i + 1);

                    if (strcmp(key, keyF) == 0)
                        F[i] = atof(value);
                    else if (strcmp(key, keyX) == 0)
                        x[i] = atof(value);
                    else if (strcmp(key, keyY) == 0)
                        y[i] = atof(value);
                    else if (strcmp(key, keyA) == 0)
                        A[i] = atof(value);
                }
            }
        }
    }
    fclose(fp);
#ifdef DEBUG_OUT_MESSAGE
    // Проверка
    printf("dx=%.3f, dy=%.3f, Nx=%d, Ny=%d\n", dx, dy, Nx, Ny);
    printf("Tx=%.2f, Ty=%.2f, область: X[%.1f; %.1f], Y[%.1f; %.1f]\n",
           Tx, Ty, Tx_min, Tx_max, Ty_min, Ty_max);
    printf("T_start=%.3f, T_step=%.3f, T_end=%.3f\n", T_start, T_step, T_end);
    printf("SNR=%s\n", SNR);
    printf("Количество сигналов: %d\n", N_Signal);

    for (int i = 0; i < N_Signal; i++)
    {
        printf("Сигнал %d: F=%.2f Hz, координаты (%.2f, %.2f), A=%.2f\n",
               i + 1, F[i], x[i], y[i], A[i]);
    }
#endif

    T = TIME_NUMBER_STEPS + 1;
    M = Nx * Ny;
    lenPhi = (Ty_max - Ty_min) / Ty + 1;
    lenTheta = (Tx_max - Tx_min) / Tx + 1;
#ifdef MULTI_THREADING
    // Инициализируем потоки
    cpu_set_t cpuset1;
    CPU_ZERO(&cpuset1);
    CPU_SET(2, &cpuset1);
    pthread_attr_setaffinity_np(&attr1, sizeof(cpu_set_t), &cpuset1);
    cpu_set_t cpuset2;
    CPU_ZERO(&cpuset2);
    CPU_SET(4, &cpuset2);
    pthread_attr_setaffinity_np(&attr2, sizeof(cpu_set_t), &cpuset2);
#endif
    InitSimulationSignal();

    InitScanScaleZone();

#ifndef OFF_ARGS
    if (strcmp(argv[1], "1D") == 0)
    {
#endif
#ifdef BEAM_1D_ON
        Init_1DArrayPattern();

        print_mem_usage();

#ifdef DEBUG_OUT_MESSAGE
        printf("Выбран режим 1D\n");
#endif

        for (int COUNTER_FOR = 0; COUNTER_FOR < NUMB_CONTER_FOR; COUNTER_FOR++)
        {
            double t = read_cpu_temp();
            if (t > 0)
                printf("Температура CPU: %.2f °C\n", t);
            else
                printf("Не удалось прочитать температуру.\n");
            start = clock();

            // computing_Rxx1D();

            compute_BEAMscan_1D();

            end = clock();
            cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG_OUT_MESSAGE
            printf("Время выполнения: %f секунд\n", cpu_time_used);
#endif
#ifdef GRAFICS

            plot_BEAMscan_x();
            plot_BEAMscan_y();
#endif
        }
#endif

#ifndef OFF_ARGS
    }
    else if (strcmp(argv[1], "2D") == 0)
    {
#endif
#ifdef BEAM_2D_ON
        Init_2DArrayPattern();

        print_mem_usage();

#ifdef DEBUG_OUT_MESSAGE
        printf("Выбран режим 2D\n");
#endif

        for (int COUNTER_FOR = 0; COUNTER_FOR < NUMB_CONTER_FOR; COUNTER_FOR++)
        {
            double t = read_cpu_temp();
            if (t > 0)
                printf("Температура CPU: %.2f °C\n", t);
            else
                printf("Не удалось прочитать температуру.\n");

            start = clock();

            computing_Rxx2D();

            compute_BEAMscan_2D();

            end = clock();
            cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
#ifdef DEBUG_OUT_MESSAGE
            printf("Время выполнения: %f секунд\n", cpu_time_used);
#endif
#ifdef GRAFICS
            plot_BEAMscan_2D();
#endif
        }
#endif
#ifndef OFF_ARGS
    }
    else
    {

        fprintf(stderr, "Ошибка: неизвестный аргумент '%s'. Используйте 1D или 2D.\n", argv[1]);

        return EXIT_FAILURE;
    }
#endif
    return EXIT_SUCCESS;
}

/**
 * @brief Инициализирует сигналы для симуляции.
 *
 * @details
 * Функция выполняет настройку необходимых обработчиков сигналов,
 * которые используются в процессе моделирования. Это может включать
 * регистрацию пользовательских функций для обработки прерываний,
 * остановки симуляции или корректного завершения программы.
 *
 * @note
 * Должна вызываться один раз в начале работы программы,
 * перед запуском основных вычислений или потоков.
 *
 * @warning
 * Неправильная инициализация сигналов может привести к некорректному
 * завершению программы или утечкам ресурсов.
 */
void InitSimulationSignal()
{
    int counter = 0;

    for (double Time = T_start; Time <= T_end; Time += T_step)
    {
        TimeSimul[counter] = Time;
        counter++;
    }

    for (int i = 0; i < N_Signal; i++)
        for (int j = 0; j <= T - 1; j++)
        {
            SinSingal[i][j].re = A[i] * cos(M_PI * 2 * F[i] * TimeSimul[j]);
            // Signal(k, tt) = (Ampl(k) * exp(1j * 2 * pi * SignalFreq(k) * t(tt)));
            SinSingal[i][j].im = A[i] * sin(M_PI * 2 * F[i] * TimeSimul[j]);
        }

#ifdef DEBUG_OUT_MATRIX_CSV

    // === Сохраняем синусоиды в файл для MATLAB ===

    FILE *f = fopen("signals.csv", "w");
    if (!f)
        return 1;

    for (int i = 0; i < N_Signal; i++)
    {
        for (int j = 0; j <= TIME_NUMBER_STEPS; j++)
        {
            fprintf(f, "%0.15f + %0.15fi",
                    SinSingal[i][j].re,
                    SinSingal[i][j].im);
            if (j < TIME_NUMBER_STEPS)
                fprintf(f, "\t"); // разделитель между отсчетами
        }
        fprintf(f, "\n"); // новая строка для следующего сигнала
    }

    fclose(f);

#endif

    for (int i = 0; i < N_Signal; i++)
    {
        int index_matrix = 0;
        for (int ny = 0; ny < Ny; ny++)
            for (int nx = 0; nx < Nx; nx++)
            {
                double angle_y = (ny - (double)(Ny - 1) / 2) * y[i] / 100.0 * M_PI;
                double angle_x = (nx - (double)(Nx - 1) / 2) * x[i] / 100.0 * M_PI;

                V[index_matrix][i].re = cos(angle_x) * cos(angle_y) - sin(angle_x) * sin(angle_y);
                V[index_matrix][i].im = sin(angle_x) * cos(angle_y) + cos(angle_x) * sin(angle_y);
                index_matrix++;
            }
    }

#ifdef DEBUG_OUT_MATRIX_CSV
    FILE *fV = fopen("V.csv", "w");
    if (!fV)
        return 1;

    for (int row = 0; row < M; row++)
    {
        for (int col = 0; col < N_Signal; col++)
        {
            fprintf(fV, "%0.15f + %0.15fi",
                    V[row][col].re,
                    V[row][col].im);
            if (col < N_Signal - 1)
                fprintf(fV, "\t");
        }
        fprintf(fV, "\n");
    }
    fclose(fV);
#endif

    // Перемножаем сигналы
    for (int t = 0; t <= T - 1; t++) // по времени
    {
        for (int m = 0; m < Nx * Ny; m++) // по элементам решётки
        {
            complex_t sum = {0, 0};
            for (int s = 0; s < N_Signal; s++) // по сигналам
            {
                // берем транспонированное: SinSignal.' → conj + поменять местами индексы
                complex_t a = SinSingal[s][t]; // Signal(s,t)
                complex_t b = V[m][s];         // V(m,s)

                sum = cadd(sum, cmul(a, b));
            }
            SignalArray[t][m] = sum;
        }
    }
}

/**
 * @brief Вычисляет одномерную автокорреляционную функцию Rxx вдоль оси X.
 *
 * @details
 * Функция реализует алгоритм вычисления автокорреляции данных
 * по пространственной координате X.
 * Используется для анализа свойств сигнала/поля и получения
 * статистических характеристик модели.
 *
 * @note
 * Предполагается, что входные данные уже загружены и подготовлены
 * для вычислений (например, заданы размеры сетки Nx, значения dx и массив данных).
 *
 * @warning
 * Некорректная инициализация параметров сетки (Nx, dx) или отсутствие
 * входных данных приведёт к неверным результатам автокорреляции.
 */
#ifdef MULTI_THREADING
#if 0
void* computing_Rxx1D_X(void* arg) {
    (void)arg;
#else
void computing_Rxx1D_X()
{
#endif
    // Основной цикл по времени
    for (int t = 0; t < T; t++)
    {
        // Вытаскиваем snapshot размером Nx × Ny
        // SignalArray.SignalArray[t][m], где m = y*Nx + x
        // (y — индекс по Y, x — индекс по X)
        // Snapshot(x,y) = SignalArray[t][y*Nx + x]

        // === Rxx_ux (по X) ===
        for (int x1 = 0; x1 < Nx; x1++)
        {
            for (int x2 = 0; x2 < Nx; x2++)
            {
                complex_t sum = {0.0, 0.0};
                for (int y = 0; y < Ny; y++)
                {
                    complex_t v1 = SignalArray[t][y * Nx + x1];
                    complex_t v2 = SignalArray[t][y * Nx + x2];
                    sum = cadd(sum, cmul(v1, cconj(v2)));
                }
                Rxx_ux[x1][x2] = cadd(Rxx_ux[x1][x2], sum);
            }
        }
    }

#ifdef DEBUG_OUT_MATRIX_CSV
    // === Запись Rxx_ux ===
    FILE *fux = fopen("Rxx_ux.csv", "w");
    if (fux)
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Nx; j++)
            {
                fprintf(fux, "%0.15f + %0.15fi",
                        Rxx_ux[i][j].re,
                        Rxx_ux[i][j].im);
                if (j < Nx - 1)
                    fprintf(fux, "\t");
            }
            fprintf(fux, "\n");
        }
        fclose(fux);
    }
#endif
}

/**
 * @brief Вычисляет одномерную автокорреляционную функцию Rxx вдоль оси Y.
 *
 * @details
 * Функция реализует расчёт автокорреляции данных по пространственной
 * координате Y. Это позволяет анализировать статистические свойства
 * моделируемого сигнала/поля в вертикальном направлении и получать
 * характеристики распределения энергии или флуктуаций.
 *
 * @note
 * Перед вызовом функции необходимо корректно инициализировать параметры
 * сетки (Ny, dy) и подготовить массив данных.
 *
 * @warning
 * Если размеры сетки заданы неправильно или отсутствуют входные данные,
 * результаты автокорреляции будут некорректными.
 */
#if 0
void* computing_Rxx1D_Y(void* arg) {
    (void)arg;
#else
void computing_Rxx1D_Y()
{
#endif
    // Основной цикл по времени
    for (int t = 0; t < T; t++)
    {
        // Вытаскиваем snapshot размером Nx × Ny
        // SignalArray.SignalArray[t][m], где m = y*Nx + x
        // (y — индекс по Y, x — индекс по X)
        // Snapshot(x,y) = SignalArray[t][y*Nx + x]

        // === Rxx_uy (по Y) ===
        for (int y1 = 0; y1 < Ny; y1++)
        {
            for (int y2 = 0; y2 < Ny; y2++)
            {
                complex_t sum = {0.0, 0.0};
                for (int x = 0; x < Nx; x++)
                {
                    complex_t v1 = SignalArray[t][y1 * Nx + x];
                    complex_t v2 = SignalArray[t][y2 * Nx + x];
                    sum = cadd(sum, cmul(cconj(v1), v2)); // поменяли местами сопряжение
                }
                Rxx_uy[y1][y2] = cadd(Rxx_uy[y1][y2], sum);
            }
        }
    }
#ifdef DEBUG_OUT_MATRIX_CSV
    // === Запись Rxx_uy ===
    FILE *fuy = fopen("Rxx_uy.csv", "w");
    if (fuy)
    {
        for (int i = 0; i < Ny; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                fprintf(fuy, "%0.15f + %0.15fi",
                        Rxx_uy[i][j].re,
                        Rxx_uy[i][j].im);
                if (j < Ny - 1)
                    fprintf(fuy, "\t");
            }
            fprintf(fuy, "\n");
        }
        fclose(fuy);
    }
#endif
}

#endif

/**
 * @brief Вычисляет одномерную автокорреляционную функцию Rxx для всей области моделирования.
 *
 * @details
 * Функция объединяет вычисления одномерных автокорреляционных функций
 * вдоль обеих осей (X и Y).
 * Обычно используется как универсальный вызов, который автоматически
 * выполняет анализ данных в обоих направлениях и сохраняет/выводит результаты.
 *
 * @note
 * Должна вызываться после инициализации параметров сетки (Nx, Ny, dx, dy)
 * и заполнения массивов исходных данных.
 *
 * @warning
 * Если входные данные не подготовлены или параметры сетки заданы неверно,
 * результаты корреляционного анализа будут недостоверны.
 *
 * @see computing_Rxx1D_X(), computing_Rxx1D_Y()
 */
void computing_Rxx1D()
{

#ifdef MULTI_THREADING


    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);

    // --- чистим атрибуты ---
    pthread_attr_destroy(&attr1);
    pthread_attr_destroy(&attr2);

#else

    // Обнуляем матрицы
    // for (int i = 0; i < Nx; i++)
    //    for (int j = 0; j < Nx; j++) {
    //        Rxx_ux[i][j].re = 0.0;
    //        Rxx_ux[i][j].im = 0.0;
    //    }

    // for (int i = 0; i < Ny; i++)
    //     for (int j = 0; j < Ny; j++) {
    //         Rxx_uy[i][j].re = 0.0;
    //         Rxx_uy[i][j].im = 0.0;
    //     }

    // Основной цикл по времени
    for (int t = 0; t < T; t++)
    {
        // Вытаскиваем snapshot размером Nx × Ny
        // SignalArray.SignalArray[t][m], где m = y*Nx + x
        // (y — индекс по Y, x — индекс по X)
        // Snapshot(x,y) = SignalArray[t][y*Nx + x]

        // === Rxx_ux (по X) ===
        for (int x1 = 0; x1 < Nx; x1++)
        {
            for (int x2 = 0; x2 < Nx; x2++)
            {
                complex_t sum = {0.0, 0.0};
                for (int y = 0; y < Ny; y++)
                {
                    complex_t v1 = SignalArray[t][y * Nx + x1];
                    complex_t v2 = SignalArray[t][y * Nx + x2];
                    sum = cadd(sum, cmul(v1, cconj(v2)));
                }
                Rxx_ux[x1][x2] = cadd(Rxx_ux[x1][x2], sum);
            }
        }

        // === Rxx_uy (по Y) ===
        for (int y1 = 0; y1 < Ny; y1++)
        {
            for (int y2 = 0; y2 < Ny; y2++)
            {
                complex_t sum = {0.0, 0.0};
                for (int x = 0; x < Nx; x++)
                {
                    complex_t v1 = SignalArray[t][y1 * Nx + x];
                    complex_t v2 = SignalArray[t][y2 * Nx + x];
                    sum = cadd(sum, cmul(cconj(v1), v2)); // поменяли местами сопряжение
                }
                Rxx_uy[y1][y2] = cadd(Rxx_uy[y1][y2], sum);
            }
        }
    }

#if 0
    // Нормировка
    double normX = (double)(T * Nx);
    double normY = (double)(T * Ny);

    for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Nx; j++) {
            Rxx_ux[i][j].re /= normX;
            Rxx_ux[i][j].im /= normX;
        }

    for (int i = 0; i < Ny; i++)
        for (int j = 0; j < Ny; j++) {
            Rxx_uy[i][j].re /= normY;
            Rxx_uy[i][j].im /= normY;
        }
#endif

#ifdef DEBUG_OUT_MATRIX_CSV
    // === Запись Rxx_ux ===
    FILE *fux = fopen("Rxx_ux.csv", "w");
    if (fux)
    {
        for (int i = 0; i < Nx; i++)
        {
            for (int j = 0; j < Nx; j++)
            {
                fprintf(fux, "%0.15f + %0.15fi",
                        Rxx_ux[i][j].re,
                        Rxx_ux[i][j].im);
                if (j < Nx - 1)
                    fprintf(fux, "\t");
            }
            fprintf(fux, "\n");
        }
        fclose(fux);
    }

    // === Запись Rxx_uy ===
    FILE *fuy = fopen("Rxx_uy.csv", "w");
    if (fuy)
    {
        for (int i = 0; i < Ny; i++)
        {
            for (int j = 0; j < Ny; j++)
            {
                fprintf(fuy, "%0.15f + %0.15fi",
                        Rxx_uy[i][j].re,
                        Rxx_uy[i][j].im);
                if (j < Ny - 1)
                    fprintf(fuy, "\t");
            }
            fprintf(fuy, "\n");
        }
        fclose(fuy);
    }
#endif
#endif
}

/**
 * @brief Вычисляет двумерную автокорреляционную функцию Rxx.
 *
 * @details
 * Функция реализует расчёт двумерной автокорреляции моделируемого поля/сигнала
 * по пространственным координатам X и Y. Такой анализ позволяет оценивать
 * пространственные зависимости и статистические свойства системы сразу в двух измерениях,
 * что важно для комплексных моделей и многомерных сигналов.
 *
 * @note
 * Перед вызовом функции необходимо корректно задать параметры сетки
 * (Nx, Ny, dx, dy) и подготовить массив данных для анализа.
 *
 * @warning
 * Некорректная инициализация сетки или отсутствие исходных данных
 * приведёт к неверным результатам вычислений.
 *
 * @see computing_Rxx1D(), computing_Rxx1D_X(), computing_Rxx1D_Y()
 */
void computing_Rxx2D()
{
    // int T = TIME_NUMBER_STEPS; // максимальный индекс времени

    // 1. Первая транспозиция: B.' (M x T+1)
    complex_t B_T[M][T + 1];
    for (int i = 0; i < M; i++)
    {
        for (int t = 0; t <= T; t++)
        {
            B_T[i][t] = SignalArray[t][i]; // строки = элементы, столбцы = время
        }
    }

    // 2. Rxx = B_T * (B_T)^H
    for (int i = 0; i < M; i++)
    { // строки Rxx
        for (int j = 0; j < M; j++)
        { // столбцы Rxx
            complex_t sum = {0.0, 0.0};
            for (int t = 0; t <= T - 1; t++)
            {
                complex_t a = B_T[i][t];     // i-я строка первой транспонированной
                complex_t b = B_T[j][t];     // j-я строка первой транспонированной
                b.im = -b.im;                // Hermitian
                sum = cadd(sum, cmul(a, b)); // аккумулируем
            }
            Rxx[i][j].re = sum.re / (double)(T);
            Rxx[i][j].im = sum.im / (double)(T);
        }
    }
#ifdef DEBUG_OUT_MATRIX_CSV

    FILE *fRxx = fopen("Rxx.csv", "w");
    if (!fRxx)
        return 1;

    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            fprintf(fRxx, "%0.15f + %0.15fi",
                    Rxx[i][j].re,
                    Rxx[i][j].im);
            if (j < M - 1)
                fprintf(fRxx, "\t");
        }
        fprintf(fRxx, "\n");
    }
    fclose(fRxx);
#endif
}

/**
 * @brief Инициализирует область сканирования и масштабирования.
 *
 * @details
 * Функция выполняет настройку параметров области моделирования,
 * в пределах которой будут выполняться вычисления или анализ сигнала.
 * Обычно используется для задания размеров зоны, шагов сетки и границ
 * по координатам X и Y. Также может подготавливать внутренние массивы
 * и структуры данных для последующего сканирования.
 *
 * @note
 * Должна вызываться перед запуском вычислительных функций
 * (например, computing_Rxx1D() или computing_Rxx2D()).
 *
 * @warning
 * Если параметры зоны заданы некорректно, вычисления в дальнейшем
 * будут давать неверные результаты или могут привести к выходу за пределы массивов.
 */
void InitScanScaleZone()
{

    int counter = 0;
    for (float i = Tx_min; i <= Tx_max; i = i + Tx)
    {
        Thetta_range[counter] = i;
        counter++;
    }
    counter = 0;
    for (float i = Ty_min; i <= Ty_max; i += Ty)
    {
        Phi_range[counter] = i;
        counter++;
    }
    // return STATE_OK;
}

/**
 * @brief Инициализирует диаграмму направленности двумерной антенной решётки.
 *
 * @details
 * Функция подготавливает двумерный массив, описывающий пространственное
 * распределение амплитудно-фазовых характеристик антенной решётки.
 * Используется для моделирования диаграммы направленности (pattern)
 * и расчёта её параметров, таких как уровень боковых лепестков, ширина
 * главного лепестка и усиление в заданных направлениях.
 *
 * @note
 * Вызывать эту функцию необходимо до выполнения процедур анализа или
 * вычисления корреляционных функций, так как диаграмма направленности
 * является исходными данными для симуляции.
 *
 * @warning
 * Некорректная инициализация может привести к ошибкам в моделировании
 * радиолокационных или коммуникационных систем, использующих данную решётку.
 */
void Init_2DArrayPattern()
{

    int lenTheta = (Tx_max - Tx_min) / Tx + 1;

    for (int i = 0; i < lenTheta; ++i)
    { // θ внешний (как в MATLAB)
        for (int j = 0; j < lenPhi; ++j)
        { // φ внутренний
            double ux = Thetta_range[i] / 100.0;
            double uy = Phi_range[j] / 100.0;

            for (int n = 0; n < Nx; ++n)
            {
                for (int m = 0; m < Ny; ++m)
                {
                    double phase_n = (n - (Nx - 1.0) / 2.0);
                    double phase_m = (m - (Ny - 1.0) / 2.0);
                    double angle = M_PI * (phase_n * ux + phase_m * uy);

                    // !!! индексация как в MATLAB: [phi, theta, m, n]
                    DigrArray[j][i][m][n].re = cos(angle);
                    DigrArray[j][i][m][n].im = sin(angle);
                }
            }
        }
    }

#ifdef DEBUG_OUT_MATRIX_CSV

    FILE *fD = fopen("DigrArray.csv", "w");
    if (!fD)
    {
        perror("Ошибка открытия файла DigrArray.csv");
        return;
    }

    for (int jp = 0; jp < lenPhi; ++jp)
    {
        for (int it = 0; it < lenTheta; ++it)
        {
            for (int ny = 0; ny < Ny; ++ny)
            {
                for (int nx = 0; nx < Nx; ++nx)
                {
                    fprintf(fD, "%0.15f + %0.15fi\n",
                            DigrArray[jp][it][ny][nx].re,
                            DigrArray[jp][it][ny][nx].im);
                }
            }
        }
    }

    fflush(fD);
    fclose(fD);

    printf("DigrArray.csv успешно записан.\n");

#endif
}

/**
 * @brief Инициализирует диаграмму направленности одномерной антенной решётки.
 *
 * @details
 * Функция формирует одномерный массив, описывающий диаграмму направленности (pattern)
 * линейной антенной решётки. Используется для моделирования характеристик излучения
 * или приёма вдоль одного пространственного направления.
 * Может включать задание амплитудно-фазового распределения, шага элементов
 * и количества излучателей.
 *
 * @note
 * Должна быть вызвана до вычисления характеристик решётки или анализа её
 * корреляционных свойств.
 *
 * @warning
 * Если параметры решётки заданы неверно (например, шаг элементов или количество
 * антенн), диаграмма направленности будет искажена, что приведёт к ошибкам
 * в дальнейших расчётах.
 */
void Init_1DArrayPattern()
{
    int lenTheta = (Tx_max - Tx_min) / Tx + 1;

    for (int i = 0; i < lenTheta; i++)
    {
        double theta = Thetta_range[i] / 100.0; // нормировка
        for (int n = 0; n < Nx; n++)
        {
            double phase_n = n - (Nx - 1) / 2.0;
            Vnm_x[n][i].re = cos(phase_n * theta * M_PI);
            Vnm_x[n][i].im = sin(phase_n * theta * M_PI);
        }
    }

    for (int j = 0; j < lenPhi; j++)
    {
        double phi = Phi_range[j] / 100.0; // нормировка
        for (int m = 0; m < Ny; m++)
        {
            double phase_m = m - (Ny - 1) / 2.0;
            Vnm_y[m][j].re = cos(phase_m * phi * M_PI);
            Vnm_y[m][j].im = sin(phase_m * phi * M_PI);
        }
    }
}

#ifdef MULTI_THREADING


/**
 * @brief Выполняет сканирование луча антенной решётки вдоль оси X.
 *
 * @details
 * Функция моделирует процесс электронного сканирования диаграммы направленности
 * двумерной или линейной антенной решётки в направлении оси X.
 * Обычно используется в многопоточной реализации: в качестве аргумента принимает
 * структуру с параметрами сканирования (например, шаг угла, количество итераций,
 * указатели на массивы данных) и выполняет вычисления в отдельном потоке.
 *
 * @param arg Указатель на структуру с параметрами сканирования (приводится к нужному типу внутри функции).
 *
 * @return void* Возвращает NULL либо указатель на результирующие данные
 * (в зависимости от реализации многопоточной логики).
 *
 * @note
 * Функция должна запускаться через `pthread_create()` или аналогичный механизм
 * многопоточности.
 *
 * @warning
 * - Перед вызовом необходимо корректно инициализировать все данные,
 *   передаваемые через `arg`.
 * - Ошибки в параметрах приведут к искажённым результатам сканирования
 *   или аварийному завершению потока.
 */
void *BeamScanX(void *arg)
{

#ifdef DEBUG_OUT_MATRIX_CSV
    printf("Available CPUs: %ld\n", sysconf(_SC_NPROCESSORS_ONLN));
    printf("[BeamScanX] running on CPU %d\n", sched_getcpu());
#endif

    computing_Rxx1D_X(NULL);

    int lenTheta = (Tx_max - Tx_min) / Tx + 1;
    double max_mag = 0.0;
    double mag = 0.0;
// === Расчет BEAMscan_x ===
#pragma GCC ivdep
    for (int i = 0; i < lenTheta; i++)
    {
        complex_t sum = {0.0, 0.0};
        for (int m = 0; m < Nx; m++)
        {
            complex_t tmp = {0.0, 0.0};
            for (int n = 0; n < Nx; n++)
            {
                tmp = cadd(tmp, cmul(cconj(Vnm_x[n][i]), cmul(Rxx_ux[n][m], Vnm_x[m][i])));
            }
            sum = cadd(sum, tmp);
        }
        BEAMscan_x[i] = sum;
        mag = sqrt(sum.re * sum.re + sum.im * sum.im);

        if (mag > max_mag)
            max_mag = mag;
    }
    // 2. Нормировка
    if (max_mag > 0.0)
    {
        for (int j = 0; j < lenTheta; j++)
        {
            BEAMscan_x[j].re /= max_mag;
            BEAMscan_x[j].im /= max_mag;
        }
    }
}


/**
 * @brief Выполняет сканирование луча антенной решётки вдоль оси Y.
 *
 * @details
 * Функция моделирует процесс электронного сканирования диаграммы направленности
 * двумерной антенной решётки в направлении оси Y.
 * Предназначена для многопоточного выполнения: принимает указатель на структуру
 * с параметрами сканирования (например, шаг угла, количество точек, буферы для
 * хранения результатов) и выполняет вычисления в отдельном потоке.
 *
 * @param arg Указатель на структуру с параметрами сканирования
 * (внутри функции приводится к соответствующему типу).
 *
 * @return void* Возвращает NULL либо указатель на структуру с результатами,
 * в зависимости от реализации логики многопоточности.
 *
 * @note
 * Обычно вызывается через `pthread_create()` для параллельного сканирования
 * по обеим осям (X и Y).
 *
 * @warning
 * - Все данные в `arg` должны быть корректно инициализированы до запуска потока.
 * - Неверные параметры или ошибки синхронизации могут привести к искажению
 *   диаграммы направленности или сбоям выполнения.
 *
 * @see BeamScanX()
 */
void *BeamScanY(void *arg)
{

#ifdef DEBUG_OUT_MATRIX_CSV
    printf("Available CPUs: %ld\n", sysconf(_SC_NPROCESSORS_ONLN));
    printf("[BeamScanX] running on CPU %d\n", sched_getcpu());
#endif

    double max_mag = 0.0;
    double mag = 0.0;
    computing_Rxx1D_Y(NULL);

// === Расчет BEAMscan_y ===
#pragma GCC ivdep
    for (int j = 0; j < lenPhi; j++)
    {
        complex_t sum = {0.0, 0.0};
        for (int p = 0; p < Ny; p++)
        {
            complex_t tmp = {0.0, 0.0};
            for (int q = 0; q < Ny; q++)
            {
                tmp = cadd(tmp, cmul(cconj(Vnm_y[q][j]), cmul(Rxx_uy[q][p], Vnm_y[p][j])));
            }
            sum = cadd(sum, tmp);
        }
        BEAMscan_y[j] = sum;

        mag = sqrt(sum.re * sum.re + sum.im * sum.im);

        if (mag > max_mag)
            max_mag = mag;
    }
    // 2. Нормировка
    if (max_mag > 0.0)
    {
        for (int j = 0; j < lenPhi; j++)
        {
            BEAMscan_y[j].re /= max_mag;
            BEAMscan_y[j].im /= max_mag;
        }
    }
}
#endif


/**
 * @brief Выполняет одномерное сканирование луча антенной решётки.
 *
 * @details
 * Функция рассчитывает диаграмму направленности (pattern) при сканировании
 * луча вдоль одного измерения (X или Y). Используется для анализа свойств
 * линейных и двумерных антенных решёток в сечении, определения положения
 * главного лепестка, ширины диаграммы и уровней боковых лепестков.
 *
 * @note
 * Перед вызовом необходимо инициализировать массив диаграммы направленности
 * (см. Init_1DArrayPattern() или Init_2DArrayPattern()) и задать параметры
 * сканирования (шаг угла, количество элементов решётки и др.).
 *
 * @warning
 * Ошибки в параметрах (например, шаге сканирования или числе элементов)
 * приведут к искажённым результатам или некорректному положению главного лепестка.
 *
 * @see compute_BEAMscan_2D(), BeamScanX(), BeamScanY()
 */
void compute_BEAMscan_1D()
{

#ifdef MULTI_THREADING
    // --- Атрибуты потока 1 ---
    // pthread_attr_init(&attr1);
    // pthread_attr_setschedpolicy(&attr1, SCHED_FIFO);
    // param1.sched_priority = sched_get_priority_max(SCHED_FIFO);
    // pthread_attr_setschedparam(&attr1, &param1);
    // pthread_attr_setinheritsched(&attr1, PTHREAD_EXPLICIT_SCHED);

    // --- Атрибуты потока 2 ---
    // pthread_attr_init(&attr2);
    // pthread_attr_setschedpolicy(&attr2, SCHED_FIFO);
    // param2.sched_priority = sched_get_priority_max(SCHED_FIFO) - 1; // чуть ниже
    // pthread_attr_setschedparam(&attr2, &param2);
    // pthread_attr_setinheritsched(&attr2, PTHREAD_EXPLICIT_SCHED);

    // --- Создание потоков уже с атрибутами ---
    pthread_create(&thread1, &attr1, BeamScanX, NULL);
    pthread_create(&thread2, &attr2, BeamScanY, NULL);

    pthread_join(thread1, NULL);
    pthread_join(thread2, NULL);

    // --- чистим атрибуты ---
    // pthread_attr_destroy(&attr1);
    // pthread_attr_destroy(&attr2);
#else

    // === Расчет BEAMscan_x ===
    for (int i = 0; i < lenTheta; i++)
    {
        complex_t sum = {0.0, 0.0};
        for (int m = 0; m < Nx; m++)
        {
            complex_t tmp = {0.0, 0.0};
            for (int n = 0; n < Nx; n++)
            {
                tmp = cadd(tmp, cmul(cconj(Vnm_x[n][i]), cmul(Rxx_ux[n][m], Vnm_x[m][i])));
            }
            sum = cadd(sum, tmp);
        }
        BEAMscan_x[i] = sum;
    }

    // === Расчет BEAMscan_y ===
    for (int j = 0; j < lenPhi; j++)
    {
        complex_t sum = {0.0, 0.0};
        for (int p = 0; p < Ny; p++)
        {
            complex_t tmp = {0.0, 0.0};
            for (int q = 0; q < Ny; q++)
            {
                tmp = cadd(tmp, cmul(cconj(Vnm_y[q][j]), cmul(Rxx_uy[q][p], Vnm_y[p][j])));
            }
            sum = cadd(sum, tmp);
        }
        BEAMscan_y[j] = sum;
    }

#endif

#ifdef DEBUG_OUT_MATRIX_CSV
    // Запись в CSV для проверки
    FILE *fx = fopen("BEAMscan_x.csv", "w");
    if (fx)
    {
        for (int i = 0; i < lenTheta; i++)
            fprintf(fx, "%0.15f + %0.15fi\n", BEAMscan_x[i].re, BEAMscan_x[i].im);
        fclose(fx);
    }

    FILE *fy = fopen("BEAMscan_y.csv", "w");
    if (fy)
    {
        for (int j = 0; j < lenPhi; j++)
            fprintf(fy, "%0.15f + %0.15fi\n", BEAMscan_y[j].re, BEAMscan_y[j].im);
        fclose(fy);
    }
#endif
}

/**
 * @brief Выполняет двумерное сканирование луча антенной решётки.
 *
 * @details
 * Функция моделирует процесс сканирования диаграммы направленности
 * двумерной антенной решётки по обоим пространственным осям (X и Y).
 * Результатом является полная двумерная диаграмма направленности,
 * отражающая распределение усиления или чувствительности в пространстве.
 * Используется для анализа характеристик антенн в радиолокационных
 * и телекоммуникационных системах.
 *
 * @note
 * Перед вызовом необходимо:
 * - инициализировать диаграмму направленности (см. Init_2DArrayPattern());
 * - корректно задать размеры решётки, шаг дискретизации и параметры сканирования.
 *
 * @warning
 * Ошибки в параметрах сетки или шагах сканирования приведут к искажению
 * диаграммы направленности или к неверной оценке боковых лепестков.
 *
 * @see compute_BEAMscan_1D(), BeamScanX(), BeamScanY()
 */
void compute_BEAMscan_2D()
{
    int lenTheta = (Tx_max - Tx_min) / Tx + 1;

    int M = Nx * Ny;

    double max_mag = 0.0; // для нормировки

    // 1. Вычисляем BEAMscan и ищем max
#pragma GCC ivdep
    for (int it = 0; it < lenTheta; ++it)
    {
        for (int jp = 0; jp < lenPhi; ++jp)
        {
            complex_t a[M];
            int idx = 0;
            for (int m = 0; m < Ny; ++m)
            {
                for (int n = 0; n < Nx; ++n)
                {
                    a[idx++] = DigrArray[jp][it][m][n]; // [φ, θ, m, n]
                }
            }

            complex_t sum = {0.0, 0.0};

            for (int i = 0; i < M; ++i)
            {
                complex_t ai_conj = cconj(a[i]);
                for (int j = 0; j < M; ++j)
                {
                    sum = cadd(sum, cmul(ai_conj, cmul(Rxx[i][j], a[j])));
                }
            }

            BEAMscan[jp][it] = sum;

            double mag = sqrt(sum.re * sum.re + sum.im * sum.im);
            if (mag > max_mag)
                max_mag = mag;
        }
    }

    // 2. Нормировка
    if (max_mag > 0.0)
    {
        for (int it = 0; it < lenTheta; ++it)
        {
            for (int jp = 0; jp < lenPhi; ++jp)
            {
                BEAMscan[jp][it].re /= max_mag;
                BEAMscan[jp][it].im /= max_mag;
            }
        }
    }

#ifdef DEBUG_OUT_MATRIX_CSV
    // 3. Запись в CSV
    FILE *fB = fopen("BEAMscan.csv", "w");
    if (!fB)
        return;

    for (int jp = 0; jp < lenPhi; ++jp)
    {
        for (int it = 0; it < lenTheta; ++it)
        {
            fprintf(fB, "%0.15f + %0.15fi",
                    BEAMscan[jp][it].re, BEAMscan[jp][it].im);
            if (it < lenTheta - 1)
                fprintf(fB, "\t");
        }
        fprintf(fB, "\n");
    }
    fclose(fB);
#endif
}

/**
 * @brief Выводит информацию об использовании оперативной памяти процессом.
 *
 * @details
 * Функция читает статистику памяти из файла `/proc/self/statm`
 * и отображает объем используемой и выделенной памяти текущим процессом.
 * Полезна для отладки и анализа производительности программы.
 *
 * @note Работает только в системах на базе Linux, где доступна файловая система `/proc`.
 *
 * @warning На других операционных системах (Windows, macOS) данный метод работать не будет.
 *
 * @see man proc (раздел о statm)
 */
void print_mem_usage()
{
    FILE *file = fopen("/proc/self/status", "r");
    if (!file)
        return;

    char line[256];
    while (fgets(line, sizeof(line), file))
    {
        if (strncmp(line, "VmRSS:", 6) == 0)
        { // Resident Set Size (реальная память)
            printf("%s", line);
        }
        if (strncmp(line, "VmSize:", 7) == 0)
        { // Виртуальная память
            printf("%s", line);
        }
    }
    fclose(file);
}

#ifdef GRAFICS

/**
 * @brief Строит график диаграммы направленности при сканировании по оси Y.
 *
 * @details
 * Функция визуализирует результаты двумерного или одномерного сканирования
 * луча антенной решётки в направлении оси Y. Обычно отображает зависимость
 * усиления (или нормированной мощности) от угла отклонения.
 * Используется для анализа характеристик диаграммы направленности:
 * ширины главного лепестка, уровня боковых лепестков, положения максимума.
 *
 * @note
 * Перед вызовом должны быть рассчитаны данные сканирования
 * (см. compute_BEAMscan_1D() или compute_BEAMscan_2D()).
 *
 * @warning
 * Если данные для оси Y не были предварительно вычислены или массив
 * результатов не инициализирован, функция построит некорректный график
 * либо вызовет ошибку доступа к памяти.
 *
 * @see plot_BEAMscan_x(), compute_BEAMscan_1D(), compute_BEAMscan_2D()
 */
void plot_BEAMscan_y()
{

    FILE *f = fopen("beam_y.dat", "w");
    for (int j = 0; j < lenPhi; j++)
    {
        double mag = sqrt(BEAMscan_y[j].re * BEAMscan_y[j].re +
                          BEAMscan_y[j].im * BEAMscan_y[j].im);
        fprintf(f, "%f %f\n", Phi_range[j], mag);
    }
    fclose(f);

    // создаем скрипт для gnuplot
    FILE *gp = popen("gnuplot -persist", "w");
    if (gp)
    {
        fprintf(gp, "set title 'BEAMscan_y'\n");
        fprintf(gp, "set xlabel 'Phi (градусы)'\n");
        fprintf(gp, "set ylabel 'Amplitude'\n");
        fprintf(gp, "plot 'beam_y.dat' with lines title 'BEAMscan_y'\n");
        fflush(gp);
        // pclose(gp); // если хочешь закрыть окно автоматически
    }
}

/**
 * @brief Строит график диаграммы направленности при сканировании по оси X.
 *
 * @details
 * Функция визуализирует результаты сканирования луча антенной решётки
 * в направлении оси X. Отображает зависимость усиления или нормированной
 * мощности от угла отклонения по горизонтальной оси.
 * Используется для анализа характеристик диаграммы направленности:
 * положения главного лепестка, ширины диаграммы, уровней боковых лепестков
 * и симметрии распределения.
 *
 * @note
 * Перед вызовом должны быть рассчитаны данные для оси X
 * (см. compute_BEAMscan_1D() или compute_BEAMscan_2D()).
 *
 * @warning
 * Если массив данных для оси X не инициализирован или повреждён,
 * построенный график будет некорректным либо произойдёт ошибка доступа к памяти.
 *
 * @see plot_BEAMscan_y(), compute_BEAMscan_1D(), compute_BEAMscan_2D()
 */
void plot_BEAMscan_x()
{

    FILE *f = fopen("beam_x.dat", "w");
    for (int j = 0; j < lenTheta; j++)
    {
        double mag = sqrt(BEAMscan_x[j].re * BEAMscan_x[j].re +
                          BEAMscan_x[j].im * BEAMscan_x[j].im);
        fprintf(f, "%f %f\n", Thetta_range[j], mag);
    }
    fclose(f);

    // создаем скрипт для gnuplot
    FILE *gp = popen("gnuplot -persist", "w");
    if (gp)
    {
        fprintf(gp, "set title 'BEAMscan_x'\n");
        fprintf(gp, "set xlabel 'Thetta (градусы)'\n");
        fprintf(gp, "set ylabel 'Amplitude'\n");
        fprintf(gp, "plot 'beam_x.dat' with lines title 'BEAMscan_x'\n");
        fflush(gp);
        // pclose(gp); // если хочешь закрыть окно автоматически
    }
}

/**
 * @brief Строит двумерный график диаграммы направленности антенной решётки.
 *
 * @details
 * Функция визуализирует результаты сканирования луча антенной решётки
 * в двухмерной плоскости (θ и φ). Отображает зависимость амплитуды
 * (или нормированной мощности) от углового положения по горизонтальной
 * и вертикальной осям. Используется для анализа характеристик
 * диаграммы направленности: положения главного лепестка, ширины диаграммы,
 * уровней боковых лепестков и симметрии распределения.
 *
 * @note
 * Перед вызовом должны быть рассчитаны данные для двумерного сканирования
 * (см. compute_BEAMscan_2D()). Функция автоматически вычисляет модуль
 * комплексных значений BEAMscan и нормирует данные.
 *
 * @warning
 * Если массив BEAMscan не инициализирован или повреждён, построенный
 * график будет некорректным либо произойдёт ошибка доступа к памяти.
 *
 * @see plot_BEAMscan_x(), plot_BEAMscan_y(), compute_BEAMscan_2D()
 */
void plot_BEAMscan_2D()
{
    FILE *f = fopen("beam_2D.dat", "w");
    if (!f)
    {
        perror("Cannot open file for 2D plot");
        return;
    }

    // Записываем данные: theta phi magnitude
    for (int it = 0; it < lenTheta; ++it)
    {
        for (int jp = 0; jp < lenPhi; ++jp)
        {
            double mag = sqrt(BEAMscan[jp][it].re * BEAMscan[jp][it].re +
                              BEAMscan[jp][it].im * BEAMscan[jp][it].im);
            fprintf(f, "%f %f %f\n", Thetta_range[it], Phi_range[jp], mag);
        }
        fprintf(f, "\n"); // пустая строка для корректного отображения сетки
    }
    fclose(f);

    FILE *gp = popen("gnuplot -persist", "w");
    if (gp)
    {
        fprintf(gp, "set title 'BEAMscan 2D'\n");
        fprintf(gp, "set xlabel 'Theta (градусы)'\n");
        fprintf(gp, "set ylabel 'Phi (градусы)'\n");
        fprintf(gp, "set zlabel 'Amplitude'\n");

        // Цветовая карта "jet" (синяя -> зеленая -> желтая -> красная)
        fprintf(gp, "set palette defined (0 'blue', 1 'cyan', 2 'green', 3 'yellow', 4 'red')\n");
        fprintf(gp, "set pm3d map\n");
        fprintf(gp, "splot 'beam_2D.dat' using 1:2:3 with pm3d notitle\n");

        fflush(gp);
        // pclose(gp);
    }
}

#endif
