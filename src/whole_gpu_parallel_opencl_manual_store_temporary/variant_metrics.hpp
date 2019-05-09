#ifndef VARIANT_METRICS_HPP
#define VARIANT_METRICS_HPP

#include <configure.hpp>
#include <omp.h>

#ifdef ENABLE_VARIANT_METRICS
#define START_TIMER(variable) (variable) = omp_get_wtime();
#define STOP_TIMER(variable) (variable) = omp_get_wtime() - (variable);
#define TIMEIT(variable, body) \
  variable = omp_get_wtime(); \
  body \
  variable = omp_get_wtime() - variable;
#else
#define TIMEIT(variable, body) body
#define START_TIMER(variable)
#define STOP_TIMER(variable)
#endif


typedef struct struct_Variant_Metrics {
  double elapsed_216;
  double elapsed_338;
  double elapsed_416;
  double elapsed_551;
  double elapsed_551_reduce;
  double elapsed_setup;
  double elapsed_compile;
  double elapsed_copy_host_to_device;
  double elapsed_copy_device_to_host;
  double elapsed_exec_setup;
  double elapsed_exec;
  double elapsed_teardown;
} Variant_Metrics;

#endif
