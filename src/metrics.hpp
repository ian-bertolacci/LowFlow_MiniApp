#ifndef METRICS_HPP
#define METRICS_HPP

#include <stdio.h>
#include <variant_metrics.hpp>
#include <util.hpp>

#ifdef ENABLE_METRICS

#define START_TIMER(variable) (variable) = omp_get_wtime();
#define STOP_TIMER(variable) (variable) = omp_get_wtime() - (variable);

#define TIMEIT( variable, ... ) \
  double MAKE_TEMP_VARIABLE_NAME_IN_LINE = omp_get_wtime(); \
  __VA_ARGS__; \
  variable = omp_get_wtime() - MAKE_TEMP_VARIABLE_NAME_IN_LINE;

#define TIMEIT_ACCUMULATE( variable, ... ) \
  double MAKE_TEMP_VARIABLE_NAME_IN_LINE = omp_get_wtime(); \
  __VA_ARGS__; \
  variable += omp_get_wtime() - MAKE_TEMP_VARIABLE_NAME_IN_LINE;

#else

#define TIMEIT( variable, ... ) __VA_ARGS__;
#define TIMEIT_ACCUMULATE( variable, ... ) __VA_ARGS__;
#define START_TIMER(variable) /* Nothing */
#define STOP_TIMER(variable) /* Nothing */

#endif


typedef struct struct_Standard_Metrics {
  double total_grid_allocation_time;
  double total_grid_deallocation_time;
  double total_grid_population_time;
} Standard_Metrics;

void printMetricInformation( FILE* stream, Standard_Metrics* standard_metrics, Variant_Metrics* metrics );
void printStandardMetricInformation( FILE* stream, Standard_Metrics* metrics );
void printVariantMetricInformation( FILE* stream, Variant_Metrics* metrics );

#endif
