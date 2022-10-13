MPIF90      ?= mpif90
INCLUDE     ?=
F90FLAGS    ?=
F90FLAGS    += -g -O2 -Wall -Wextra -pedantic -fbounds-check -Wuninitialized \
               -Wno-unused-parameter -cpp -ffree-form  -ffree-line-length-none $(USEBLAS)

# OPENBLAS stuff
OPENBLASBASE    ?= #/public/apps/openblas/gcc-4.8.5
BLASLIBS        ?= -fexternal-blas -L$(OPENBLAS)/lib -lopenblas
F90FLAGS        += -DOPENBLAS $(BLASLIBS)

# HDF5
HDF5     ?= #/public/apps/HDF5/1.13.2_gcc-4.8.5
HDF5INC  += -I$(HDF5)/include
F90FLAGS += $(HDF5INC) -L$(HDF5)/lib -lhdf5_fortran -lhdf5
