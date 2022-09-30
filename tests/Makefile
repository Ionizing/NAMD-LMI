include ../make.inc

MOD = ../src
LIB = $(MOD)/libnamd.a
BIN = test.x

test: $(BIN)
	for b in $(BIN); do ./$${b}; done

fruit.o: fruit.f90

%.x: %.f90 $(LIB) fruit.o
	$(MPIF90) -I$(MOD) $< -o $@

clean:
	rm -f *.x *.o *.mod
