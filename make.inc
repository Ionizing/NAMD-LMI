MPIF90      ?= mpif90
INCLUDE     ?=
F90FLAGS    ?=
F90FLAGS    += -g -O2 -fPIC -Wall -Wextra -pedantic -fbounds-check -Wuninitialized \
               -Wno-unused-parameter -cpp -ffree-form  -ffree-line-length-none \
               -frecord-gcc-switches

# OPENBLAS stuff
OPENBLAS    ?= #/public/apps/openblas/gcc-4.8.5
BLASLIBS    ?= -fexternal-blas -L$(OPENBLAS)/lib -lopenblas
F90FLAGS    += -DOPENBLAS $(BLASLIBS)

# MKL stuff
#MKLFLAGS  = -fexternal-blas -m64 -I"${MKLROOT}/include" -L${MKLROOT}/lib/intel64 -lmkl_rt -Wl,--no-as-needed -lpthread -lm -ldl
#F90FLAGS += $(MKLFLAGS)

# HDF5
HDF5     ?= #/public/apps/HDF5/1.13.2_gcc-4.8.5
HDF5INC  += -I$(HDF5)/include
F90FLAGS += $(HDF5INC) -DHDF5 -L$(HDF5)/lib -lhdf5_fortran -lhdf5

GIT_VERSION := "$(shell git describe --dirty --always --tags)"
F90FLAGS += -DGIT_VERSION=\"$(GIT_VERSION)\"
