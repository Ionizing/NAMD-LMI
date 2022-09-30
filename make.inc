.SUFFIXES: .o .f90

MPIF90      ?= mpif90
F90FLAGS    ?= -g -O2 -Wall -Wextra -pedantic -fbounds-check -Wuninitialized \
               -Wno-unused-parameter -cpp

.f90.o:
	$(MPIF90) $(F90FLAGS) -c $<
