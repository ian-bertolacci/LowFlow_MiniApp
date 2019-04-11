cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# List of variants (names are also directories)
set(VARIANTS simple_whole_serial simple_whole_cpu_parallel_omp)
# List all source files
set(BASE_SOURCE_FILES global_types.cpp main.cpp)

# Create rooted path
PREPEND(BASE_SOURCE_FILES ${BASE_SOURCE_DIR} ${BASE_SOURCE_FILES})

# Add local include
include_directories(${BASE_SOURCE_DIR})

#Enable debug mode globally
set(CMAKE_BUILD_TYPE Debug)

# Execute all variant cmake files
foreach(VARIANT ${VARIANTS})
  add_subdirectory( ${VARIANT} )
endforeach(VARIANT)
