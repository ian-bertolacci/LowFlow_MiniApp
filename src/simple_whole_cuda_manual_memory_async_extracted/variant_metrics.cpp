#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  fprintf( stream,
    "Elapsed temp_alloc:          %f\n"
    "Elapsed temp_dealloc:        %f\n",
    //"Elapsed setup:               %f\n"
    //"Elapsed copy_host_to_device: %f\n"
    //"Elapsed copy_device_to_host: %f\n"
    //"Elapsed exec_setup:          %f\n"
    //"Elapsed exec:                %f\n"
    //"Elapsed teardown:            %f\n",
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
