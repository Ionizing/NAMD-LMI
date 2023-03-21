#define TINT2STR(x) TRIM(int2str(x))
#define TFINT2STR(x,fmt) TRIM(int2str(x,fmt))
#define TREAL2STR(x) TRIM(real2str(x))
#define AT "at " // __FILE__ // ":" // TINT2STR(__LINE__)

#ifdef __INTEL_COMPILER
#define REALPART DREAL
#define IMAGPART DIMAG
#else
#define STOP CALL BACKTRACE; stop
#endif


!! Profiling routines
#ifdef PROFILING
#define PROFILING_DECLARATIONS      INTEGER(8) :: profiling_start, profiling_rate, profiling_end
#define PROFILING_START             CALL SYSTEM_CLOCK(profiling_start)
#define PROFILING_END(msg)          CALL SYSTEM_CLOCK(profiling_end, profiling_rate); \
                                    PRINT '("[PROFILING of ", A,", @", A, I0, ":] ", F12.4, " (ms)")', \
                                        msg, __FILE__, __LINE__, \
                                        DBLE(profiling_end - profiling_start) * 1E3 / profiling_rate
#else
#define PROFILING_DECLARATIONS
#define PROFILING_START
#define PROFILING_END(msg)
#endif
