include ../make.inc

OBJS = string.o common.o wavecar.o
LIB = libnamd.a

namd_lumi.x: namd_lumi.f90 $(LIB)
	$(MPIF90) -o $@ $< $(F90FLAGS)

$(LIB): $(OBJS) lib.o
	ar -rc $@ $^

lib.o: $(OBJS)

wavecar.o: wavecar.f90 common.h common.o

clean:
	rm -f *.x *.o *.mod *.a

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $<
