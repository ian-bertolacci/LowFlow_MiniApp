#ifndef VARIANT_CONFIGURE_HPP
#define VARIANT_CONFIGURE_HPP

#include <global_configure.hpp>
#include <Kokkos_Core.hpp>

typedef struct struct_VariantOptions {

} VariantOptions;

static const VariantOptions VariantOptions_Default;

static const struct option variant_long_options[] = {
  {"help", no_argument, nullptr, 'h'},
  LONG_OPTIONS_TERMINAL
};

/*
Possible options for kokkos_execution_space:
+ USE_SERIAL: Use Kokkos::Serial
+ USE_OPENMP:  Use Kokkos::OpenMP
+ USE_PTHREADS:  Use Kokkos::PThreads
+ USE_CUDA:  Use Kokkos::CUDA
*/

// Check that at most one is defined
#if defined(USE_SERIAL) + defined(USE_OPENMP) + defined(USE_PTHREADS) + defined(USE_CUDA) <= 1
  // Define correct execution space
  #if defined(USE_SERIAL)
    typedef Kokkos::Serial kokkos_execution_space;
    typedef Kokkos::HostSpace kokkos_storage_space;
  #elif defined(USE_OPENMP)
    typedef Kokkos::OpenMP kokkos_execution_space;
    typedef Kokkos::HostSpace kokkos_storage_space;
  #elif defined(USE_PTHREADS)
    typedef Kokkos::Threads kokkos_execution_space;
    typedef Kokkos::HostSpace kokkos_storage_space;
  #elif defined(USE_CUDA)
    typedef Kokkos::Cuda kokkos_execution_space;
    typedef Kokkos::HostSpace kokkos_host_storage_space;
    typedef Kokkos::CudaSpace kokkos_device_storage_space;
    typedef kokkos_device_storage_space kokkos_storage_space;
  #else
    typedef Kokkos::Serial kokkos_execution_space;
    typedef Kokkos::HostSpace kokkos_storage_space;
  #endif
// Error out if more than one is defined
#else
  #error "At most one of the USE_<Execution Space> Kokkos options may be may be defined!"
  // Report which ones are used.
  #if defined(USE_SERIAL)
    #error "USE_SERIAL is defined"
  #endif
  #if defined(USE_OPENMP)
    #error "USE_OPENMP is defined"
  #endif
  #if defined(USE_PTHREADS)
    #error "USE_PTHREADS is defined"
  #endif
  #if defined(USE_CUDA)
    #error "USE_CUDA is defined"
  #endif
#endif

#endif
