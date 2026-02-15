# Makefile — сборка в текущей директории
CC ?= $(CROSS_COMPILE)gcc
export CC
ARCH_FLAGS ?= -march=native

# main по-умолчанию
MAIN ?= main.c
TARGET = $(basename $(MAIN))

SRC := $(MAIN) $(COMMON_SRC)
OBJ := $(patsubst %.c,%.o,$(SRC))

# --- Опции сборки (y/n) ---
DEBUG_OUT_MESSAGE ?= y
DEBUG_OUT_MATRIX_CSV ?= n
MULTI_THREADING ?= n
BEAM_2D_ON ?= y
BEAM_1D_ON ?= y
GRAFICS ?= n

# --- Формирование DEFS (используем := для немедленного вычисления) ---
DEFS := $(if $(filter y,$(DEBUG_OUT_MESSAGE)),-DDEBUG_OUT_MESSAGE,) \
        $(if $(filter y,$(DEBUG_OUT_MATRIX_CSV)),-DDEBUG_OUT_MATRIX_CSV,) \
        $(if $(filter y,$(MULTI_THREADING)),-DMULTI_THREADING,) \
        $(if $(filter y,$(BEAM_2D_ON)),-DBEAM_2D_ON,) \
        $(if $(filter y,$(BEAM_1D_ON)),-DBEAM_1D_ON,) \
        $(if $(filter y,$(GRAFICS)),-DGRAFICS,) \
        -DTHETTA_NUMBER_POINT=241 \
        -DPHI_NUMBER_POINT=241 \
        -DNUMBER_X=5 \
        -DNUMBER_Y=6 \
        -DNUMBER_SIGNALS=2 \
        -DSNAPSHOTS=5 \
        -DNUMB_CONTER_FOR=100 \
        -DMAX_LINE=256

# --- Флаги компиляции с DEFS ---
CFLAGS = -Wall -Wextra -Ofast $(ARCH_FLAGS) -ffast-math -funroll-loops -fdiagnostics-color=always $(DEFS)
LDFLAGS = -lm -lpthread

.PHONY: all clean check-defs

all: $(TARGET)

# --- Компиляция объектов с выводом команды для отладки ---
%.o: %.c
	@echo "Компиляция $< с DEFS: $(DEFS)"
	$(CC) $(CFLAGS) -c $< -o $@

# --- Линковка ---
$(TARGET): $(OBJ)
	$(CC) $(OBJ) -o $@ $(LDFLAGS)

# --- Очистка ---
clean:
	rm -f $(OBJ) $(TARGET)
	rm -f *.csv *.dat

# --- Проверка DEFS ---
check-defs:
	@echo "DEFS = $(DEFS)"
	@echo "CFLAGS = $(CFLAGS)"
	@echo ""
	@echo "Первые 20 символов DEFS:"
	@echo "$(DEFS)" | cut -c1-20
	@echo ""
	@echo "Команда компиляции для main.c:"
	@echo "$(CC) $(CFLAGS) -c main.c -o main.o"