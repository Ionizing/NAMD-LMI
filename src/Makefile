include ../make.inc

SRCS = $(filter-out namd_lumi.f90,$(wildcard *.f90))
OBJS = $(filter-out lib.o,$(SRCS:.f90=.o))
LIB = libnamd.a
LIBEXTERNAL = ../external/libexternal.a

namd_lumi.x: namd_lumi.f90 $(LIB) $(LIBEXTERNAL)
	$(MPIF90) -o $@ $^ $(F90FLAGS) -I../external

$(LIB): $(OBJS) lib.o
	ar -rcv $@ $^

lib.o: $(OBJS)

wavecar.o: string.o wavecar.f90 common.h common.o

tdm.o: string.o common.o wavecar.o

nac.o: string.o common.o wavecar.o input.o tdm.o

hamiltonian.o: string.o nac.o wavecar.o common.o

input.o: string.o

surface_hopping.o: common.o string.o hamiltonian.o input.o nac.o

clean:
	rm -f *.x *.o *.mod *.a

%.o: %.f90
	$(MPIF90) -c $< $(F90FLAGS) -I../external

$(LIBEXTERNAL):
	$(MAKE) -C ../external
