# Makefile — сборка в текущей директории
CC ?= $(CROSS_COMPILE)gcc
export CC
ARCH_FLAGS ?= -march=native
CFLAGS = -Wall -Wextra -Ofast $(ARCH_FLAGS) -ffast-math -funroll-loops -fdiagnostics-color=always

LDFLAGS = -lm -lpthread

# main по-умолчанию
MAIN ?= main.c
TARGET = $(basename $(MAIN))

SRC := $(MAIN) $(COMMON_SRC)
OBJ := $(patsubst %.c,%.o,$(SRC))

# --- Опции сборки (y/n)
# Параметр вывода отладочных сообщений в консоль
DEBUG_OUT_MESSAGE ?= y
# Параметр сохранения .csv для отладки и сравнения с matlab рабочим кодом
DEBUG_OUT_MATRIX_CSV ?= n
# разделяем расчет 2x1D Beam Scan на 2 потока
MULTI_THREADING ?= y
MULTI_THREADING ?= n
# Поддержка расчета 2D Beam Scan
BEAM_2D_ON ?= y
# Поддержка расчета 1D Beam Scan
BEAM_1D_ON ?= y
# Поддержка вывода графика
GRAFICS ?= n

# --- DEFS передаются компилятору
DEFS = $(if $(filter y,$(DEBUG_OUT_MESSAGE)),-DDEBUG_OUT_MESSAGE,) \
       $(if $(filter y,$(DEBUG_OUT_MATRIX_CSV)),-DDEBUG_OUT_MATRIX_CSV,) \
       $(if $(filter y,$(MULTI_THREADING)),-DMULTI_THREADING,) \
       $(if $(filter y,$(BEAM_2D_ON)),-DBEAM_2D_ON,) \
       $(if $(filter y,$(BEAM_1D_ON)),-DBEAM_1D_ON,) \
       $(if $(filter y,$(GRAFICS)),-DGRAFICS,) \
       -DTHETTA_NUMBER_POINT=241 \
       -DPHI_NUMBER_POINT=241 \
       -DNUMBER_X=2 \
       -DNUMBER_Y=2 \
       -DNUMBER_SIGNALS=2 \
       -DSNAPSHOTS=5 \
       -DNUMB_CONTER_FOR=100 \
       -DMAX_LINE=256 \

.PHONY: all clean

# --- По умолчанию собираем бинарник
all: $(TARGET)

# --- Компиляция объектов
%.o: %.c
	$(CC) $(CFLAGS) $(DEFS) -c $< -o $@

# --- Линковка
$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

# --- Очистка
clean:
	rm -f $(OBJ) $(TARGET)
	rm -f *.csv *.dat
