cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# List of variants (names are also directories)
set(VARIANTS
  simple_whole_serial
  simple_whole_serial_agi
  simple_whole_serial_agi_access_using_domain
  simple_whole_serial_manual_fuse
  simple_whole_cpu_parallel_omp
  simple_whole_cpu_parallel_omp_manual_fuse
  whole_gpu_parallel_opencl_manual_store_temporary
  whole_gpu_parallel_opencl_manual_fuse_store_temporary
  simple_whole_kokkos
  simple_whole_cuda
  simple_whole_cuda_lambda
  simple_whole_cuda_manual_memory
  simple_whole_cuda_single_kernel
)

# List all source files
set(BASE_SOURCE_FILES global_types.cpp global_configure.cpp metrics.cpp main.cpp)

# Create rooted path
PREPENDPATH(BASE_SOURCE_FILES ${BASE_SOURCE_DIR} ${BASE_SOURCE_FILES})

# Add local include
include_directories(${BASE_SOURCE_DIR})

# Enable Release mode globally
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Compiler and compiler flags
# Set defaults first so variants can rewrite as necessary
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

if(${WITH_CUDA})
  if(${WITH_KOKKOS})
    set(CMAKE_CUDA_COMPILER nvcc_wrapper)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED")
    set(CMAKE_CUDA_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_EXTENSIONS OFF) # use std=c++11 instead of std=gnu++11
  endif()

  commafy_flags("${CMAKE_CXX_FLAGS}" __CMAKE_CXX_FLAGS_comma_separated)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda -Xcompiler ${__CMAKE_CXX_FLAGS_comma_separated}")
  set(CMAKE_CUDA_STANDARD 11)
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  if(${WITH_CUDA})
    if("${OpenMP_C_FLAGS}")
      commafy_flags("${OpenMP_C_FLAGS}" __OpenMP_C_FLAGS_comma_separated)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${__OpenMP_C_FLAGS_comma_separated}")
    endif()

    if("${OpenMP_CXX_FLAGS}")
      commafy_flags("${OpenMP_CXX_FLAGS}" __OpenMP_CXX_FLAGS_comma_separated)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler ${__OpenMP_CXX_FLAGS_comma_separated}")
    endif()

    if("${OpenMP_EXE_LINKER_FLAGS}")
      commafy_flags("${OpenMP_EXE_LINKER_FLAGS}" __OpenMP_EXE_LINKER_FLAGS_comma_separated)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xlinker ${__OpenMP_EXE_LINKER_FLAGS_comma_separated}")
    endif()
  endif()

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} -fopenmp")
else()
  message( WARNING "OpenMP is mandatory, but not found. manually setting flags." )
  if(${WITH_CUDA})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler -fopenmp -Xlinker -fopenmp")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
  endif()
endif()

# Execute all variant cmake files
foreach(VARIANT ${VARIANTS})
  # set(VARIANT_DIR ${VARIANT})
  add_subdirectory( ${VARIANT} )
endforeach(VARIANT)
