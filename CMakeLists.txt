cmake_minimum_required(VERSION 3.5)

project(LowFlow)

set(CMAKE_DIR ${PROJECT_SOURCE_DIR}/cmake)

include(${CMAKE_DIR}/utilities.txt)
include(${CMAKE_DIR}/build_variant.txt)

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})
set(INSTALL_BIN_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(BASE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)

add_subdirectory( ${BASE_SOURCE_DIR} )
