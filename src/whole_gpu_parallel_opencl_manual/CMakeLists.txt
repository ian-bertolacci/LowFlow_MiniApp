cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# Variant name
set(VARIANT_NAME simple_whole_gpu_parallel_opencl)
# List source files
set(VARIANT_SOURCE_FILES variant_types.cpp variant_science.cpp)
# List libraries
set(VARIANT_LINK_LIBRARIES OpenCL)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

CREATE_VARIANT_TARGET(${VARIANT} "${BASE_SOURCE_FILES}" "${VARIANT_SOURCE_FILES}" ${INSTALL_BIN_DIR} "${VARIANT_LINK_LIBRARIES}")
