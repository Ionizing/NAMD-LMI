#OPENBLAS    ?= /public/apps/openblas/gcc-4.8.5
#USEBLAS     ?= -fexternal-blas -L$(OPENBLAS)/lib -lopenblas
MPIF90      ?= mpif90
F90FLAGS    ?= -g -O2 -Wall -Wextra -pedantic -fbounds-check -Wuninitialized \
               -Wno-unused-parameter -cpp -ffree-form  -ffree-line-length-none $(USEBLAS)
