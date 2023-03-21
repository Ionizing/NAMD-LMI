include ../make.inc

MOD = -I../src -I../external
INC = ../src/common.h
LIB = ../src/libnamd.a
LIBMOD = ../src/namd_lumi_mod.mod
LIBEXTERNAL = ../external/libexternal.a
BIN = test.x

TEST_SRCS = $(wildcard test_*.f90)
TEST_OBJS = $(TEST_SRCS:.f90=.o)

test: test.x
	./$<
	
test.x: test.f90 $(TEST_OBJS) $(LIB) $(LIBEXTERNAL)
	$(MPIF90) $(MOD) $^ -o $@ $(F90FLAGS)

%.o: %.f90 $(LIBEXTERNAL) $(INC) $(LIBMOD)
	$(MPIF90) -c $< $(MOD) $(F90FLAGS)

$(LIB) $(LIBMOD):
	$(MAKE) -C ../src

$(LIBEXTERNAL):
	$(MAKE) -C ../external

clean:
	rm -f test *.o *.mod *.x *.h5
