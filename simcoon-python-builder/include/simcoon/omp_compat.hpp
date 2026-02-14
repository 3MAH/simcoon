#pragma once

#ifdef _OPENMP
#include <omp.h>
inline void omp_set_active_levels(int levels) {
#if _OPENMP >= 201811
    omp_set_max_active_levels(levels);
#else
    omp_set_nested(levels > 0);
#endif
}
#endif
