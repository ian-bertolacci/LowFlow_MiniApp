#ifndef VARIANT_CONFIGURE_HPP
#define VARIANT_CONFIGURE_HPP

#include <global_configure.hpp>

// Variant specific metrics:
// INTEL_OPENCL : Using Intel's OpenCL runtime (or something like that)
// NVIDIA_OPENCL : Using NVidia's OpenCL runtime (or something like that)

typedef struct struct_VariantOptions {
  bool verbose;
} VariantOptions;

static const VariantOptions VariantOptions_Default {
  .verbose = false
};

static const struct option variant_long_options[] = {
  {"help",           no_argument, nullptr, 'h'},
  {"verbose",        no_argument, nullptr, 'v'},
  LONG_OPTIONS_TERMINAL
};

enum VariantExitCode : ExitType {
  OPENCL_DEVICE_CREATION_FAILURE  = 100,
  OPENCL_COMPILATION_FAILURE      = 101,
  OPENCL_KERNEL_CREATION_FAILURE  = 102,
  OPENCL_BUFFER_WRITE_FAILURE     = 103,
  OPENCL_BUFFER_READ_FAILURE      = 104,
  OPENCL_KERNEL_SETUP_FAILURE     = 105,
  OPENCL_KERNEL_EXECUTION_FAILURE = 106
};

#endif
