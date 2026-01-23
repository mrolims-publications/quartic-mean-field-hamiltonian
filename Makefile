# -------- Compiler and default flags --------
CC ?= gcc
CFLAGS = -Wall -Wextra -Wpedantic -O3

# -------- Project layout --------
PROGRAMS := programs
BINARIES := bin
SRC_FILES := $(wildcard src/*.c)

# -------- Compiler-specific flags --------
ifeq ($(CC),icx)
	CFLAGS += -ipo
else
	CFLAGS += -flto
endif

# -------- System-specific linking --------
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	# macOS (Homebrew)
	GSL_INC := -I/opt/homebrew/include
	GSL_LIB := -L/opt/homebrew/lib
	LFLAGS  := $(GSL_LIB) -lgsl -lgslcblas -lm -framework Accelerate
else
	# Linux
	GSL_INC := -I/usr/include
	GSL_LIB := -L/usr/lib
	LFLAGS  := $(GSL_LIB) -lgsl -lgslcblas -lm -llapack
endif

CFLAGS += $(GSL_INC)

ALL_BINS := \
	$(BINARIES)/std_on_site_energy_history.x \
	$(BINARIES)/lle_history.x \
	$(BINARIES)/q_history.x \
	$(BINARIES)/p_PDF.x \
	$(BINARIES)/numerical_moments.x

.PHONY: all clean \
	std_on_site_energy_history lle_history q_history p_PDF numerical_moments

all: $(ALL_BINS)

$(BINARIES):
	mkdir -p $@

$(BINARIES)/%.x: $(PROGRAMS)/%.c $(SRC_FILES) | $(BINARIES)
	$(CC) $(CFLAGS) $(SRC_FILES) $< -o $@ $(LFLAGS)

std_on_site_energy_history: $(BINARIES)/std_on_site_energy_history.x
lle_history:               $(BINARIES)/lle_history.x
q_history:                 $(BINARIES)/q_history.x
p_PDF:            $(BINARIES)/p_PDF.x
numerical_moments:         $(BINARIES)/numerical_moments.x

clean:
	rm -rf $(BINARIES)/*.x
