#pragma once

#include <cstdio>
#include <cstdlib>

#define asserta(cond)                                                          \
    do {                                                                       \
        if (!(cond)) {                                                         \
            std::fprintf(stderr, "asserta failed: %s:%d: %s\n", __FILE__,     \
                         __LINE__, #cond);                                     \
            std::abort();                                                      \
        }                                                                      \
    } while (0)
