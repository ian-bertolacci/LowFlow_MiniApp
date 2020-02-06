#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  fprintf( stream,
    "Elapsed fused:               %f\n"
    "Elapsed 551_reduce:          %f\n"
    "Elapsed temp_grid_alloc:     %f\n"
    "Elapsed temp_grid_dealloc:   %f\n",
    //"Elapsed setup:               %f\n"
    //"Elapsed copy_host_to_device: %f\n"
    //"Elapsed copy_device_to_host: %f\n"
    //"Elapsed exec_setup:          %f\n"
    //"Elapsed exec:                %f\n"
    //"Elapsed teardown:            %f\n",
    metrics->elapsed_fused,
    metrics->elapsed_551_reduce,
    metrics->elapsed_temp_grid_alloc,
    metrics->elapsed_temp_grid_dealloc
    //metrics->elapsed_setup,
    //metrics->elapsed_copy_host_to_device,
    //metrics->elapsed_copy_device_to_host,
    //metrics->elapsed_exec_setup,
    //metrics->elapsed_exec,
    //metrics->elapsed_teardown
  );
}
