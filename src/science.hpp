#ifndef SCIENCE_HPP
#define SCIENCE_HPP
#include <types.hpp>
#include <configure.hpp>
#include <metrics.hpp>

void science(
  Variant_Domain* domain,    // Iteration Domain

  Variant_Grid* fp,          // Output grid fp
  Variant_Grid* vx,          // Output grid vx
  Variant_Grid* vy,          // Output grid vy
  Variant_Grid* vz,          // Output grid vz

  Variant_Grid* dp,          // Input grid dp
  Variant_Grid* et,          // Input grid et
  Variant_Grid* odp,         // Input grid odp
  Variant_Grid* opp,         // Input grid opp
  Variant_Grid* osp,         // Input grid osp
  Variant_Grid* permxp,      // Input grid permxp
  Variant_Grid* permyp,      // Input grid permyp
  Variant_Grid* permzp,      // Input grid permzp
  Variant_Grid* pop,         // Input grid pop
  Variant_Grid* pp,          // Input grid pp
  Variant_Grid* rpp,         // Input grid rpp
  Variant_Grid* sp,          // Input grid sp
  Variant_Grid* ss,          // Input grid ss
  Variant_Grid* z_mult_dat,  // Input grid z_mult_dat
  Variant_Grid* x_ssl_dat,   // Input grid x_ssl_dat
  Variant_Grid* y_ssl_dat,   // Input grid y_ssl_dat

  VariantOptions options,    // Variant specific runtime configuration options
  Variant_Metrics* metrics   // Output metrics for variant
);

#endif
