#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  fprintf( stream,
    "Elapsed 216:                 %fs\n"
    "Elapsed 338:                 %fs\n"
    "Elapsed 416:                 %fs\n"
    "Elapsed 551:                 %fs\n"
    "Elapsed 551_reduce:          %fs\n"
    "Elapsed setup:               %fs\n"
    "Elapsed copy_host_to_device: %fs\n"
    "Elapsed copy_device_to_host: %fs\n"
    "Elapsed exec_setup:          %fs\n"
    "Elapsed exec:                %fs\n"
    "Elapsed teardown:            %fs\n",
    metrics->elapsed_216,
    metrics->elapsed_338,
    metrics->elapsed_416,
    metrics->elapsed_551,
    metrics->elapsed_551_reduce,
    metrics->elapsed_setup,
    metrics->elapsed_copy_host_to_device,
    metrics->elapsed_copy_device_to_host,
    metrics->elapsed_exec_setup,
    metrics->elapsed_exec,
    metrics->elapsed_teardown
  );
}
