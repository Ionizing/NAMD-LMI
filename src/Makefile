include ../make.inc

OBJS = lib.o
LIB = libnamd.a

namd_lumi.x: namd_lumi.f90 $(LIB)
	$(MPIF90) -o $@ $< $(F90FLAGS)

libnamd.a: $(OBJS)
	ar -rc $@ $<

clean:
	rm -f *.x *.o *.mod *.a
