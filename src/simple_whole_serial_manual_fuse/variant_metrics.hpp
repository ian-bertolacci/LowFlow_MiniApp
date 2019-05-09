#ifndef VARIANT_METRICS_HPP
#define VARIANT_METRICS_HPP

#include <configure.hpp>
#include <omp.h>

#ifdef ENABLE_VARIANT_METRICS
#define TIMEIT(variable, body) \
  variable = omp_get_wtime(); \
  body \
  variable = omp_get_wtime() - variable;
#else
#define TIMEIT(variable, body) body
#endif

typedef struct struct_Variant_Metrics {
  double elapsed_fused;
  double elapsed_551;
} Variant_Metrics;

#endif
