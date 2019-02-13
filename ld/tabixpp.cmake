cmake_minimum_required(VERSION 3.5)
project(tabixpp VERSION 0.0.1)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)

if(CGET_PREFIX)
  include_directories(${CGET_PREFIX}/include)
endif()

if(CGET_PREFIX)
  link_directories(${CGET_PREFIX}/lib)
endif()

set(MAKEFILE_COMPILER_FLAGS "-g -Wall -O2 -fPIC -D_FILE_OFFSET_BITS=64 -D_USE_KNETFILE")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MAKEFILE_COMPILER_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MAKEFILE_COMPILER_FLAGS}")

find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
add_library(tabixpp
  tabix.hpp
  tabix.cpp)
target_link_libraries(tabixpp ${HTS_LIB})
install(TARGETS tabixpp DESTINATION lib)
install(FILES tabix.hpp DESTINATION include RENAME tabixpp.hpp)
