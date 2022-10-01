include ../make.inc

OBJS = string.o
LIB = libnamd.a

namd_lumi.x: namd_lumi.f90 $(LIB)
	$(MPIF90) -o $@ $< $(F90FLAGS)

$(LIB): $(OBJS) lib.o
	ar -rc $@ $^

lib.o: $(OBJS)

clean:
	rm -f *.x *.o *.mod *.a
