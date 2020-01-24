#include <metrics.hpp>

void printMetricInformation( FILE* stream, Standard_Metrics* standard_metrics, Variant_Metrics* variant_metrics ){
  #ifdef ENABLE_METRICS
    printStandardMetricInformation( stream, standard_metrics );
    printVariantMetricInformation( stream, variant_metrics );
  #else
    /* Do nothing */
    fprintf( stream, "Metrics disabled.\n" );
  #endif
}

void printStandardMetricInformation( FILE* stream, Standard_Metrics* metrics ){
  fprintf( stream,
    "Total grid allocation time:   %f\n"
    "Total grid deallocation time: %f\n"
    "Total grid population time:   %f\n",
    metrics->total_grid_allocation_time,
    metrics->total_grid_deallocation_time,
    metrics->total_grid_population_time
  );
}
