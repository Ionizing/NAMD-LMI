include ../make.inc

SOURCE = $(wildcard *.f90)
OBJS   = $(patsubst %.f90,%.o,$(SOURCE))
LIB    = libexternal.a

$(LIB): $(OBJS)
	ar -rcv $@ $^

clean:
	rm -f *.o *.mod *.a

fitpack.o: fitpack_core.o

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $<
