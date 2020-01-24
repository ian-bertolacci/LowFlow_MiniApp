cmake_minimum_required(VERSION 3.5)
include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

# If using CUDA
if(${WITH_CUDA})
  set(VARIANT_NAME simple_whole_cuda)
  set(VARIANT_SOURCE_FILES variant_types.cpp variant_science.cpp variant_metrics.cpp variant_configure.cpp variant_setup_teardown.cpp)

  # an annoying consequence of CMake matching the compiler to the file extension
  # is that cpp files, which may even include some CUDA code, are compiled by
  # the C++ compiler, and not the CUDA compiler. This is fixed by setting these
  # files as being in the CUDA language. May also set them to be in the CUDA
  # language for other variants; however, that has not yet been seen in practise.
  set_source_files_properties(${VARIANT_SOURCE_FILES} PROPERTIES LANGUAGE CUDA)
  set_source_files_properties(${BASE_SOURCE_FILES} PROPERTIES LANGUAGE CUDA)

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
    CONFIG_MACRO_DEFNS ENABLE_METRICS
  )
endif()
