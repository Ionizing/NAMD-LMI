include ../make.inc

OBJS = string.o common.o wavecar.o tdm.o nac.o hamiltonian.o input.o
LIB = libnamd.a
LIBEXTERNAL = ../external/libexternal.a

namd_lumi.x: namd_lumi.f90 $(LIB) $(LIBEXTERNAL)
	$(MPIF90) -o $@ $^ $(F90FLAGS) -I../external

$(LIB): $(OBJS) lib.o
	ar -rcv $@ $^

lib.o: $(OBJS)

wavecar.o: string.o wavecar.f90 common.h common.o

tdm.o: string.o common.o wavecar.o

nac.o: string.o common.o wavecar.o

hamiltonian.o: string.o nac.o wavecar.o common.o

input.o: string.o

clean:
	rm -f *.x *.o *.mod *.a

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $<

$(LIBEXTERNAL):
	$(MAKE) -C ../external
