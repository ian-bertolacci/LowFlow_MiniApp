cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# Variant name
set(VARIANT_NAME simple_whole_cpu_parallel_omp_manual_fuse)
# List source files
set(VARIANT_SOURCE_FILES variant_types.cpp variant_science.cpp variant_metrics.cpp variant_configure.cpp variant_setup_teardown.cpp)
# List libraries
set(VARIANT_LINK_LIBRARIES)

CREATE_VARIANT_TARGET(${VARIANT} "${BASE_SOURCE_FILES}" "${VARIANT_SOURCE_FILES}" ${INSTALL_BIN_DIR} "${VARIANT_LINK_LIBRARIES}" "")
CREATE_VARIANT_TARGET(${VARIANT} "${BASE_SOURCE_FILES}" "${VARIANT_SOURCE_FILES}" ${INSTALL_BIN_DIR} "${VARIANT_LINK_LIBRARIES}" "ENABLE_VARIANT_METRICS")
