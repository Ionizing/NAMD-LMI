#define TINT2STR(x) TRIM(int2str(x))
#define TREAL2STR(x) TRIM(real2str(x))
#define AT "at " // __FILE__ // ":" // TINT2STR(__LINE__)
#define STOP CALL BACKTRACE; stop
