include ../make.inc

MOD = -I../src -I../external
INC = ../src/common.h
LIB = ../src/libnamd.a
LIBEXTERNAL = ../external/libexternal.a
BIN = test.x

.PHONY: $(LIB)

TEST_OBJS = $(patsubst %.f90,%.o,$(wildcard test_*.f90))

test: test.f90 $(TEST_OBJS) $(LIB) $(LIBEXTERNAL)
	echo $(TEST_OBJS)
	$(MPIF90) $(MOD) $^ -o $@.x $(F90FLAGS)
	./$@.x

%.o: %.f90 $(LIBEXTERNAL) $(INC) ../src/namd_lumi_mod.mod
	$(MPIF90) -c $< $(MOD) $(F90FLAGS)

$(LIB):
	$(MAKE) -C ../src

$(LIBEXTERNAL):
	$(MAKE) -C ../external

clean:
	rm -f test *.o *.mod *.x *.h5
