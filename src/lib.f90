MODULE namd_lumi_mod
    USE string_mod
    USE common_mod
    USE input_mod
    USE wavecar_mod
    USE tdm_mod
    USE nac_mod
    USE hamiltonian_mod
    USE surface_hopping_mod

    !< Bridge to scripts.c
    !< Provides:
    !<     void print_scripts_gen_efield()
    !<     void print_scripts_plot()
    INTERFACE
        SUBROUTINE print_scripts_gen_efield() BIND(C, NAME='print_scripts_gen_efield')
        END SUBROUTINE print_scripts_gen_efield
        SUBROUTINE print_scripts_plot()       BIND(C, NAME='print_scripts_plot')
        END SUBROUTINE print_scripts_plot
    END INTERFACE
END MODULE namd_lumi_mod
