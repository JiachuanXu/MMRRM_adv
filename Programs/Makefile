# C compiler and flags
FFTW3LIB  = -lfftw3
MATHLIB   = -lm
GSLLIB    = -lgsl -lgslcblas
GSLINC    = -I/usr/local/include
FFTW3INC  = -I/usr/local/include
CFLAGS    = -I/usr/local/include
CC        = gcc -Wall -g -fopenmp

# Directory containing header files
HEADER_DIR= /header

# Object files
OBJ_FILES = test_powers test_mudecomp test_pdf test_statool test_power_ratio

# Non-project dependence
HEADERS   = ${HEADER_DIR}/mmrrm_adv.h \
	${HEADER_DIR}/read_args.h \
	${HEADER_DIR}/read_write_box.h

test_MMRRM: MMRRM_adv.c \
	${OBJ_FILES} \

	${CC} ${CFLAGS} ${FFTW3INC} -o test_MMRRM MMRRM_adv.c ${FFTW3LIB} ${MATHLIB}

test_powers: powers.c 

	${CC} ${CFLAGS} ${FFTW3INC} -o test_powers powers.c ${FFTW3LIB} ${MATHLIB}

test_statool: statool.c

	${CC} ${CFLAGS} ${FFTW3INC} -o test_statool statool.c ${FFTW3LIB} ${MATHLIB}

test_pdf: pixelfluxpdf.c

	${CC} ${CFLAGS} ${FFTW3INC} -o test_pdf pixelfluxpdf.c ${FFTW3LIB} ${MATHLIB}

test_mudecomp: mudecomp.c

	${CC} ${CFLAGS} ${FFTW3INC} ${GSLINC} -o test_mudecomp mudecomp.c ${FFTW3LIB} ${GSLLIB} ${MATHLIB}

test_power_ratio: power_ratio.c

	${CC} ${CFLAGS} ${FFTW3INC} -o test_power_ratio power_ratio.c ${FFTW3LIB} ${MATHLIB}

#test_quasi_linear: quasi_linear.c

#	${CC} ${CFLAGS} ${FFTW3INC} -o test_quasi_linear quasi_linear.c ${FFTW3LIB} ${MATHLIB}

.PHONY:clean
clean:
	rm test_MMRRM ${OBJ_FILES}
