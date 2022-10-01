include ../make.inc

MOD = ../src
INC = $(MOD)/common.h
LIB = $(MOD)/libnamd.a
BIN = test.x

.PHONY: $(LIB)

TEST_OBJS = test_lib.o
FRUIT = fruit.o

test: test.f90 $(TEST_OBJS) $(LIB) $(FRUIT)
	$(MPIF90) -I$(MOD) $^ -o $@ $(F90FLAGS)
	./$@

%.o: %.f90 $(FRUIT) $(INC)
	$(MPIF90) -c $< -I$(MOD) $(F90FLAGS)

$(LIB):
	$(MAKE) -C ../src

clean:
	rm -f test *.o *.mod
