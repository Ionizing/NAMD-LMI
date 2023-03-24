#include <stdio.h>

#define STR2(x) #x
#define STR(x) STR2(x)

#ifdef __APPLE__
#define USTR(x) "_" STR(x)
#else
#define USTR(x) STR(x)
#endif

#ifdef _WIN32
#define INCBIN_SECTION ".rdata, \"dr\""
#elif defined __APPLE__
#define INCBIN_SECTION "__TEXT,__const"
#else
#define INCBIN_SECTION ".rodata"
#endif

// this aligns start address to 16 and terminates byte array with explict 0
// which is not really needed, feel free to change it to whatever you want/need
#define INCBIN(prefix, name, file) \
    __asm__(".section " INCBIN_SECTION "\n" \
            ".global " USTR(prefix) "_" STR(name) "_start\n" \
            ".balign 16\n" \
            USTR(prefix) "_" STR(name) "_start:\n" \
            ".incbin \"" file "\"\n" \
            \
            ".global " STR(prefix) "_" STR(name) "_end\n" \
            ".balign 1\n" \
            USTR(prefix) "_" STR(name) "_end:\n" \
            ".byte 0\n" \
    ); \
    extern __attribute__((aligned(16))) const char prefix ## _ ## name ## _start[]; \
    extern                              const char prefix ## _ ## name ## _end[];

INCBIN(script, gen_efield, "./scripts/gen_efield.py")
INCBIN(script, plot,       "./scripts/plot.py")


/*
 * In fact I can write a macro to do the repetitive work, but for the sake of simplicity,
 * I choose the more verbose way. However if more scripts need to be written in the 
 * following way, writing a macro to do the work is possible.
 */


void print_scripts_gen_efield() {
    const char fname[] = "gen_efield.py";
    const int len      =   (char*)&script_gen_efield_end
                         - (char*)&script_gen_efield_start;

    fprintf(stdout, "[INFO] Writting pre-processing script to %s ...\n", fname);

    FILE* f = fopen(fname, "w");
    if (NULL == f) {
        fprintf(stderr, "[ERROR] Open file %s failed!", fname);
    }

    fwrite(script_gen_efield_start, sizeof (char), len, f);
    fprintf(stdout, "[INFO] %s written.\n", fname);
    fclose(f);
}


void print_scripts_plot() {
    const char fname[] = "plot.py";
    const int len      =   (char*)&script_plot_end
                         - (char*)&script_plot_start;

    fprintf(stdout, "[INFO] Writting post-processing script to %s ...\n", fname);

    FILE* f = fopen(fname, "w");
    if (NULL == f) {
        fprintf(stderr, "[ERROR] Open file %s failed!", fname);
    }

    fwrite(script_plot_start, sizeof (char), len, f);
    fprintf(stdout, "[INFO] %s written.\n", fname);
    fclose(f);
}
