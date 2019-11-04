#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
#ifdef ENABLE_VARIANT_METRICS
  fprintf( stream,
    "Elapsed 216:        %fs\n"
    "Elapsed 338:        %fs\n"
    "Elapsed 416:        %fs\n"
    "Elapsed 551:        %fs\n"
    "Elapsed 551_reduce: %fs\n",
    metrics->elapsed_216,
    metrics->elapsed_338,
    metrics->elapsed_416,
    metrics->elapsed_551,
    metrics->elapsed_551_reduce
  );
#else
  /* Do nothing */
  fprintf( stream, "Variant metric disabled.\n" );
#endif
}
