#ifndef VARIANT_CONFIGURE_HPP
#define VARIANT_CONFIGURE_HPP

#include <global_configure.hpp>

typedef struct struct_VariantOptions {
  int chunks;
  int streams;
} VariantOptions;

static const VariantOptions VariantOptions_Default {
  .chunks = 10,
  .streams = 5
};

static const struct option variant_long_options[] = {
  {"chunks",  required_argument, nullptr, 'c'},
  {"streams", required_argument, nullptr, 'r'},
  {"help",    no_argument,       nullptr, 'h'},
  LONG_OPTIONS_TERMINAL
};

#endif
