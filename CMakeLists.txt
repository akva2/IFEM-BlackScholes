project(BlackScholes)

cmake_minimum_required(VERSION 2.6)

# Add local modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${PROJECT_SOURCE_DIR}/../../cmake/Modules)

if(NOT IFEM_CONFIGURED)
  find_package(IFEM REQUIRED)
  include_directories(${IFEM_INCLUDE_DIRS})
  add_definitions(${IFEM_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${IFEM_CXX_FLAGS}")
endif()

set(BS_SRCS BlackScholes.C)

include_directories(${IFEM_INCLUDES} ${PROJECT_SOURCE_DIR} ../Common)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_library(CommonBS STATIC ${BS_SRCS})

add_executable(BlackScholes main_BlackScholes.C)
target_link_libraries(BlackScholes CommonBS IFEMAppCommon ${IFEM_LIBRARIES})
list(APPEND CHECK_SOURCES ${BS_SRCS} main_BlackScholes.C)

# Installation
install(TARGETS BlackScholes DESTINATION bin COMPONENT bin)

# For generating the doxy
set(EXTRA_DOXY_PATHS "${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR}")
#add_doc_target(BlackScholes BlackScholes)

# Regression tests
enable_testing()
include(IFEMTesting)
if(MPI_FOUND)
  # Add parallel tests here
endif()
if(NOT MPI_FOUND OR IFEM_SERIAL_TESTS_IN_PARALLEL)
endif()
list(APPEND TEST_APPS BlackScholes)

if(IFEM_COMMON_APP_BUILD)
  set(TEST_APPS ${TEST_APPS} PARENT_SCOPE)
  set(UNIT_TEST_NUMBER ${UNIT_TEST_NUMBER} PARENT_SCOPE)
else()
  add_check_target()
endif()
