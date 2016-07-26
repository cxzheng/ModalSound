#ifndef CONFIG_H
#   define CONFIG_H

#ifdef USE_DOUBLE_PRECI
    typedef double REAL;
#elif USE_SINGLE_PRECI
    typedef float REAL;
#else
    typedef double REAL;
#endif

#endif
