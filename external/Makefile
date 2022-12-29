include ../make.inc

SOURCE = $(wildcard *.f90)
OBJS   = $(patsubst %.f90,%.o,$(SOURCE))
LIB    = libexternal.a

$(LIB): $(OBJS)
	ar -rcv $@ $^

clean:
	rm -f *.o *.mod *.a

fftpack5.o: fftpack5.f90
	$(MPIF90) $(filter-out -pedantic,$(F90FLAGS)) -c $<

%.o: %.f90
	$(MPIF90) $(F90FLAGS) -c $<
