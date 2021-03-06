#ifndef VARIANT_METRICS_HPP
#define VARIANT_METRICS_HPP

#include <configure.hpp>
#include <omp.h>

typedef struct struct_Variant_Metrics {
  double elapsed_216;
  double elapsed_338;
  double elapsed_416;
  double elapsed_551;
  double elapsed_551_reduce;
  double elapsed_setup;
  double elapsed_copy_host_to_device;
  double elapsed_copy_device_to_host;
  double elapsed_exec_setup;
  double elapsed_exec;
  double elapsed_teardown;
} Variant_Metrics;

#endif
