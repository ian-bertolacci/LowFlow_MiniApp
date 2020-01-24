#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
#ifdef ENABLE_METRICS
  fprintf( stream,
    "Elapsed 216:        %f\n"
    "Elapsed 338:        %f\n"
    "Elapsed 416:        %f\n"
    "Elapsed 551:        %f\n"
    "Elapsed 551_reduce: %f\n",
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
