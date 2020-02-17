#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  fprintf( stream,
    "Elapsed 216 prepare:         %f\n"
    "Elapsed 216:                 %f\n"
    "Elapsed 338 prepare:         %f\n"
    "Elapsed 338:                 %f\n"
    "Elapsed 416 prepare:         %f\n"
    "Elapsed 416:                 %f\n"
    "Elapsed 551 prepare:         %f\n"
    "Elapsed 551:                 %f\n"
    "Elapsed 551_reduce:          %f\n"
    "Elapsed copyback:            %f\n"
    "Elapsed free device:         %f\n"
    "Elapsed temp_alloc:          %f\n"
    "Elapsed temp_dealloc:        %f\n",
    //"Elapsed setup:               %f\n"
    //"Elapsed copy_host_to_device: %f\n"
    //"Elapsed copy_device_to_host: %f\n"
    //"Elapsed exec_setup:          %f\n"
    //"Elapsed exec:                %f\n"
    //"Elapsed teardown:            %f\n",
    metrics->elapsed_prepare_216,
    metrics->elapsed_216,
    metrics->elapsed_prepare_338,
    metrics->elapsed_338,
    metrics->elapsed_prepare_416,
    metrics->elapsed_416,
    metrics->elapsed_prepare_551,
    metrics->elapsed_551,
    metrics->elapsed_551_reduce,
    metrics->elapsed_copyback,
    metrics->elapsed_free_device,
    metrics->elapsed_temp_alloc,
    metrics->elapsed_temp_dealloc
    //metrics->elapsed_setup,
    //metrics->elapsed_copy_host_to_device,
    //metrics->elapsed_copy_device_to_host,
    //metrics->elapsed_exec_setup,
    //metrics->elapsed_exec,
    //metrics->elapsed_teardown
  );
}
