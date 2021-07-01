cmake_minimum_required(VERSION 3.9)
set(CMAKE_VERBOSE_MAKEFILE ON)

include(ExternalProject)

if(CGET_PREFIX)
  include_directories(${CGET_PREFIX}/include)
endif()

if(CGET_PREFIX)
  link_directories(${CGET_PREFIX}/lib)
endif()

message(STATUS "Current build directory: " ${CMAKE_CURRENT_BINARY_DIR})
message(STATUS "Current source directory: " ${CMAKE_CURRENT_SOURCE_DIR})
message(STATUS "CMake install libdir: " ${CMAKE_INSTALL_LIBDIR})
message(STATUS "CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "cget install prefix: " ${CGET_INSTALL_PREFIX})
message(STATUS "cget prefix: " ${CGET_PREFIX})
message(STATUS "Current build output root directory: " ${BUILD_OUTPUT_ROOT_DIRECTORY})
message(STATUS "CMake build type: " ${CMAKE_BUILD_TYPE})
message(STATUS "CMake install include dir: " ${CMAKE_INSTALL_INCLUDEDIR})

set(ARROW_PREFIX "${BUILD_OUTPUT_ROOT_DIRECTORY}")
set(ARROW_INCLUDE_DIR "${ARROW_PREFIX}/include")
set(ARROW_LIB_DIR "${ARROW_PREFIX}")
set(ARROW_SHARED_LIB "${ARROW_LIB_DIR}/libarrow${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(ARROW_STATIC_LIB "${ARROW_LIB_DIR}/libarrow.a")

set(ARROW_CMAKE_ARGS
  # Build settings
  #-DARROW_BUILD_STATIC=ON
  #-DARROW_BUILD_SHARED=OFF
  #-DARROW_BOOST_USE_SHARED=ON
  -DARROW_BUILD_TESTS=OFF
  -DARROW_OPTIONAL_INSTALL=ON
  #-DARROW_TEST_MEMCHECK=OFF
  #-DARROW_BUILD_BENCHMARKS=OFF

  # Arrow modules/dependencies
  #-DARROW_WITH_LZ4=ON
  -DARROW_WITH_ZSTD=ON
  #-DARROW_WITH_BROTLI=ON
  #-DARROW_WITH_SNAPPY=ON
  #-DARROW_WITH_ZLIB=ON
  #-DARROW_FLIGHT=ON
  #-DARROW_HIVESERVER2=ON
  #-DARROW_ORC=ON
  #-DARROW_GANDIVA=ON
  #-DARROW_GANDIVA_JAVA=ON
  -DARROW_PARQUET=ON
  -DARROW_FILESYSTEM=ON
  #-DARROW_HDFS=ON
  #-DARROW_IPC=ON
  ##-DARROW_COMPUTE=OFF
  #-DARROW_CUDA=OFF
  #-DARROW_GPU=OFF
  ##-DARROW_JEMALLOC=OFF
  ##-DARROW_BOOST_VENDORED=OFF
  #-DARROW_PYTHON=ON
)

add_custom_target(arrow
  ALL
  COMMAND cd cpp && mkdir -p ${CMAKE_BUILD_TYPE} && cd ${CMAKE_BUILD_TYPE} && ${CMAKE_COMMAND} -DCMAKE_TOOLCHAIN_FILE=${CGET_PREFIX}/cget/cget.cmake -DCMAKE_INSTALL_PREFIX=${CGET_PREFIX} ${ARROW_CMAKE_ARGS} .. && ${CMAKE_COMMAND} --build .
	VERBATIM
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(arrow_install
  COMMAND cd cpp/${CMAKE_BUILD_TYPE} && ${CMAKE_MAKE_PROGRAM} install
  VERBATIM
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_dependencies(arrow_install arrow)

install(CODE "
execute_process(
  COMMAND ${CMAKE_COMMAND} --build . --target arrow_install
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
")
