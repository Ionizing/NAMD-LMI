#include "../src/common.h"

MODULE test_input
    USE fruit
    USE string_mod
    USE input_mod

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_input_fn
        CALL test_input_nml
    END SUBROUTINE test_input_fn


    SUBROUTINE test_input_nml
        TYPE(input) :: inp, inp_read

        inp%rundir      = "../../run"
        inp%wavetype    = "ncl"
        inp%ikpoint     = 2
        inp%brange      = [100, 200]
        inp%nsw         = 3000
        inp%ndigit      = 5
        inp%dt          = 0.5
        inp%nsample     = 5
        inp%inispin     = 1

        ALLOCATE(inp%inibands(inp%nsample))
        ALLOCATE(inp%inisteps(inp%nsample))

        inp%inibands = [100, 101, 103, 105, 200]
        inp%inisteps = [  1, 200, 500, 100, 114]

        CALL input_to_file(inp, fname="test_input.nml")
        !CALL input_to_file(inp)
        CALL input_from_file(inp_read, fname="test_input.nml")

        CALL assert_equals(inp%rundir,    inp_read%rundir,    AT)
        CALL assert_equals(inp%wavetype,  inp_read%wavetype,  AT)
        CALL assert_equals(inp%ikpoint,   inp_read%ikpoint,   AT)
        CALL assert_equals(inp%brange,    inp_read%brange, 2, AT)
        CALL assert_equals(inp%nsw,       inp_read%nsw,       AT)
        CALL assert_equals(inp%ndigit,    inp_read%ndigit,    AT)
        CALL assert_equals(inp%dt,        inp_read%dt,        AT)
        CALL assert_equals(inp%nsample,   inp_read%nsample,   AT)
        CALL assert_equals(inp%inispin,   inp_read%inispin,   AT)
        CALL assert_equals(inp%inibands,  inp_read%inibands, inp%nsample, AT)
        CALL assert_equals(inp%inisteps,  inp_read%inisteps, inp%nsample, AT)
    END SUBROUTINE test_input_nml
END MODULE test_input
