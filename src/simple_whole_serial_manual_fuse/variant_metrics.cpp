#include <metrics.hpp>

void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics ){
  fprintf( stream,
    "Elapsed fused: %f\n"
    "Elapsed 551: %f\n",
    metrics->elapsed_fused,
    metrics->elapsed_551
  );
}
