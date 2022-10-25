include ../make.inc

OBJS = string.o common.o wavecar.o tdm.o nac.o hamiltonian.o
LIB = libnamd.a

namd_lumi.x: namd_lumi.f90 $(LIB)
	$(MPIF90) -o $@ $^ $(F90FLAGS)

$(LIB): $(OBJS) lib.o
	ar -rcv $@ $^

lib.o: $(OBJS)

wavecar.o: string.o wavecar.f90 common.h common.o

tdm.o: string.o common.o wavecar.o

nac.o: string.o common.o wavecar.o

hamiltonian.o: string.o nac.o wavecar.o common.o

clean:
	rm -f *.x *.o *.mod *.a

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $<
