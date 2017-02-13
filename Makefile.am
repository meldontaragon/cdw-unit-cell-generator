SUFFIXES = .c .h

WARN_FLAGS = -Wall -Wextra -pedantic $(WARN_EXTRA)
L_FLAGS = -lm

AM_CFLAGS = $(WARN_FLAGS) -std=c99
LDADD = $(L_FLAGS)

bin_PROGRAMS = bin/generateTMDSuperCell
bin_generateTMDSuperCell_SOURCES = src/generateTMDSuperCell.c src/tmd_functions.c src/printing.c\
	src/fractional.c src/atomic_sites.c src/structure.c\
	src/testing.c src/elements.h src/header.h


debug:
	make clean; make "CFLAGS=$(AM_CFLAGS) $(ERROR_FLAG) $(O_FLAG_0)";

optim:
	make clean; make "CFLAGS=$(AM_CFLAGS) $(O_FLAG_3)";

.PHONY : debug optim

#extra flags

O_FLAG_0 = -g3 -pg -O0 
O_FLAG_1 = -g1 -pg -O1 
O_FLAG_2 = -O2
O_FLAG_3 = -O3

WARN_EXTRA = \
-Wfloat-equal\
-Wunused\
-Wunknown-pragmas\
-Wwrite-strings\
-Winit-self\
-Wcast-align\
-Wcast-qual\
-Wpointer-arith\
-Wstrict-aliasing\
-Wmissing-declarations\
-Wmissing-include-dirs\
-Wno-unused-parameter\
-Wuninitialized\
-Wcast-align\
-Wcast-qual\
-Wdisabled-optimization\
-Wformat=2\
-Wredundant-decls\
-Wshadow\
-Wstrict-overflow=4\
-Wswitch-default\
-Wundef\
-Wno-unused\
-Wsign-conversion\
-Wmissing-format-attribute\
-Wtype-limits\
-Wconversion\
-Wsign-compare\
-Wmissing-field-initializers\
-Wno-multichar\
-Winline\
$(GCC_WARN_EXTRA)

GCC_WARN_EXTRA=\
-Wsuggest-attribute=pure\
-Wsuggest-attribute=const\
-Wsuggest-attribute=noreturn\
-Wlogical-op

GCC_NEW_WARN_EXTRA=\
-Wmaybe-uninitialized\
-Wsuggest-attribute=format\
-Wnarrowing

ERROR_FLAG=\
-Werror\
-Wno-error=unused-parameter\
-Wno-error=unused-variable\
-Wno-error=pedantic

