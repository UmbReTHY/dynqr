cmake_minimum_required(VERSION 3.0)
project(dynqr VERSION 0.1 LANGUAGES CXX)

add_library(dynqr INTERFACE)
target_include_directories(dynqr INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src)

option(BUILD_TESTS "Build unit tests" OFF)
if(BUILD_TESTS)
  enable_testing()
  include(catch.cmake)

  add_executable(dynqr_test_exe dynqr_test.cc)
  target_link_libraries(dynqr_test_exe PRIVATE dynqr catch)
endif(BUILD_TESTS)