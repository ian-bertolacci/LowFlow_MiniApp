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
)

# List all source files
set(BASE_SOURCE_FILES global_types.cpp global_configure.cpp main.cpp)

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

# Flag options
if(${CMAKE_CXX_COMPILER} MATCHES "nvcc")
  set(CMAKE_CXX_FLAGS "--expt-extended-lambda -D_FORCE_INLINES -D_MWAITXINTRIN_H_INCLUDED -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
  set(CMAKE_CXX_EXTENSIONS OFF) # use std=c++11 instead of std=gnu++11
else()
  set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-unused-parameter -Wno-unused-variable")
  set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
  set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
  warning( "OpenMP is mandatory, but not found. manually setting flags." )
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# Execute all variant cmake files
foreach(VARIANT ${VARIANTS})
  # set(VARIANT_DIR ${VARIANT})
  add_subdirectory( ${VARIANT} )
endforeach(VARIANT)
