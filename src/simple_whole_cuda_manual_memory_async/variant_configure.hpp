#ifndef VARIANT_CONFIGURE_HPP
#define VARIANT_CONFIGURE_HPP

#include <global_configure.hpp>

typedef struct struct_VariantOptions {

} VariantOptions;

static const VariantOptions VariantOptions_Default;

static const struct option variant_long_options[] = {
  {"help", no_argument, nullptr, 'h'},
  LONG_OPTIONS_TERMINAL
};

#endif
