#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#include ".vscode/vscode_macros.h"

// Переменные для определения времени расчета
clock_t start, end;
double cpu_time_used;

//  В будующем может, когда будет из этого делаться модуль или драйвер для расчета
typedef enum
{
    UNKNOWN_ERROR = 0,
    OVER_CALCULATION = 1,
    PARAMS_IS_MISSING = 2,
    START = 4,
} StateMath;

// Структура для комплексного числа

typedef struct
{
    double re;
    double im;
} complex_t;


// Инициализация Массивов для расчета

complex_t SinSingal[NUMBER_SIGNALS][SNAPSHOTS];

complex_t V[NUMBER_X * NUMBER_Y][NUMBER_SIGNALS];

complex_t SignalArray[NUMBER_SIGNALS][NUMBER_X * NUMBER_Y];
complex_t Rxx_ux[NUMBER_X][NUMBER_X];
complex_t Rxx_uy[NUMBER_Y][NUMBER_Y];
complex_t Vnm_x[NUMBER_X][THETTA_NUMBER_POINT];
complex_t Vnm_y[NUMBER_Y][PHI_NUMBER_POINT];
double BEAMscan_x[THETTA_NUMBER_POINT]; // одномерный спектр по theta
double BEAMscan_y[PHI_NUMBER_POINT];    // одномерный спектр по phi

complex_t DigrArray[PHI_NUMBER_POINT][THETTA_NUMBER_POINT][NUMBER_Y][NUMBER_X];
complex_t Rxx[NUMBER_X * NUMBER_Y][NUMBER_X * NUMBER_Y];
double BEAMscan[PHI_NUMBER_POINT][THETTA_NUMBER_POINT]; // сохраняем комплекс

float Thetta_range[THETTA_NUMBER_POINT];
float Phi_range[PHI_NUMBER_POINT];   

// переменные для потокового распределения

#ifdef MULTI_THREADING
    pthread_t thread1, thread2;
    pthread_attr_t attr1, attr2;
    struct sched_param param1, param2;
#endif

// Основные функции
void InitSimulationSignal();
void computing_Rxx1D();
void computing_Rxx2D();
void InitScanScaleZone();
void Init_2DArrayPattern();
void Init_1DArrayPattern();
void compute_BEAMscan_2D();
void compute_BEAMscan_1D();
void plot_BEAMscan_y();
void plot_BEAMscan_x();
void plot_BEAMscan_2D();
void print_mem_usage();