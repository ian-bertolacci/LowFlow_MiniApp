cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

set(VARIANT_NAME simple_whole_cpu_parallel_omp_manual_fuse)
set(VARIANT_SOURCE_FILES variant_types.cpp variant_science.cpp variant_metrics.cpp variant_configure.cpp variant_setup_teardown.cpp)

CREATE_VARIANT_TARGET(
  NAME ${VARIANT}
  BASE_SRC_FILES ${BASE_SOURCE_FILES}
  VARIANT_SRC_FILES ${VARIANT_SOURCE_FILES}
  INSTALL_DIR ${INSTALL_BIN_DIR}
)

CREATE_VARIANT_TARGET(
  NAME ${VARIANT}
  BASE_SRC_FILES ${BASE_SOURCE_FILES}
  VARIANT_SRC_FILES ${VARIANT_SOURCE_FILES}
  INSTALL_DIR ${INSTALL_BIN_DIR}
  CONFIG_MACRO_DEFNS ENABLE_VARIANT_METRICS
)
