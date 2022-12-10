MPIF90      ?= mpif90
INCLUDE     ?=
F90FLAGS    ?=
F90FLAGS    += -g -O2 -fPIC -Wall -Wextra -pedantic -fbounds-check -Wuninitialized \
               -Wno-unused-parameter -cpp -ffree-form  -ffree-line-length-none \
               -frecord-gcc-switches -funroll-all-loops

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

GIT_HASH    := $(shell git describe --match=nEvErMaTcH --always --abbrev=8 --dirty)
GIT_VERSION := $(shell git describe --dirty --always --tags | sed -nE 's/^v([0-9]+)(.([0-9]+)(.([0-9]+))?)?(.*)$$/\1 \3 \5/p')

VER_MAJOR   := $(or $(word 1, $(GIT_VERSION)), 0)
VER_MINOR   := $(or $(word 2, $(GIT_VERSION)), 0)
VER_PATCH   := $(or $(word 3, $(GIT_VERSION)), 0)

F90FLAGS += -DGIT_HASH=\"$(GIT_HASH)\" \
            -DVER_MAJOR=$(VER_MAJOR) \
            -DVER_MINOR=$(VER_MINOR) \
            -DVER_PATCH=$(VER_PATCH)
