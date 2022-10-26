include ../make.inc

MOD = ../src
INC = $(MOD)/common.h
LIB = $(MOD)/libnamd.a
BIN = test.x

.PHONY: $(LIB)

TEST_OBJS = $(patsubst %.f90,%.o,$(wildcard test_*.f90))
FRUIT = fruit.o

test: test.f90 $(TEST_OBJS) $(LIB) $(FRUIT)
	$(MPIF90) -I$(MOD) $^ -o $@.x $(F90FLAGS)
	./$@.x

%.o: %.f90 $(FRUIT) $(INC) $(LIB)
	$(MPIF90) -c $< -I$(MOD) $(F90FLAGS)

$(LIB):
	$(MAKE) -C ../src

clean:
	rm -f test *.o *.mod *.x *.h5
