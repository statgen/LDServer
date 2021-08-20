find_package(benchmark REQUIRED)
find_library(BOOST_SYSTEM_LIB boost_system HINT ${CGET_PREFIX}/lib)
find_library(BOOST_IOSTREAMS_LIB boost_iostreams HINT ${CGET_PREFIX}/lib)
find_library(MKL_RT_LIB mkl_rt)
find_package(Threads REQUIRED)

find_package(Python3 3.8 EXACT COMPONENTS Development REQUIRED)

find_library(OPENBLAS_LIB openblas HINT ${CGET_PREFIX}/lib)
message(STATUS "OpenBLAS = ${OPENBLAS_LIB}")

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

if(NOT MKL_RT_LIB)
   set(MKL_RT_LIB "")
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(benchmark-ldserver
  benchmark.cpp)

target_link_libraries(benchmark-ldserver
  LDServer
  benchmark::benchmark
  ${Python3_LIBRARIES} 
  ${BOOST_SYSTEM_LIB}
  ${BOOST_IOSTREAMS_LIB}
  ${OPENBLAS_LIB}
  ${MKL_RT_LIB}
  ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS benchmark-ldserver COMPONENT cli RUNTIME DESTINATION test OPTIONAL)
