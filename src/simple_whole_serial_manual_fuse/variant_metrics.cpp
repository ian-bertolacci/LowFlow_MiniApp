#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
#ifdef ENABLE_VARIANT_METRICS
  fprintf( stream,
    "Elapsed fused: %f\n"
    "Elapsed 551: %f\n",
    metrics->elapsed_fused,
    metrics->elapsed_551
  );
#else
  /* Do nothing */
  fprintf( stream, "Variant metric disabled.\n" );
#endif
}
