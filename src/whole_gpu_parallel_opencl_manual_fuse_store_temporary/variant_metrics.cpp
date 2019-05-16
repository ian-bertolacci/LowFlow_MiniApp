#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  #ifdef ENABLE_VARIANT_METRICS
    fprintf( stream,
      "Elapsed 216 338 416 551 fuse: %f\n"
      "Elapsed 551 reduce: %f\n"
      "Elapsed setup: %f\n"
      "Elapsed compile: %f\n"
      "Elapsed host->device: %f\n"
      "Elapsed device->host: %f\n"
      "Elapsed setup execution: %f\n"
      "Elapsed execution: %f\n"
      "Elapsed teardown: %f\n",

      metrics->elapsed_216_338_416_551,
      metrics->elapsed_551_reduce,
      metrics->elapsed_setup,
      metrics->elapsed_compile,
      metrics->elapsed_copy_host_to_device,
      metrics->elapsed_copy_device_to_host,
      metrics->elapsed_exec_setup,
      metrics->elapsed_exec,
      metrics->elapsed_teardown
    );
  #else
    /* Do nothing */
    fprintf( stream, "Variant metric disabled.\n" );
  #endif
}
