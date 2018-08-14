cmake_minimum_required(VERSION 3.8)
project(LDServer VERSION 0.0.1)
include(GNUInstallDirs)
add_compile_options(-Wall -Wpedantic)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

find_library(SAVVY_LIB savvy HINT ${CGET_PREFIX}/lib)
find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(ZSTD_LIB zstd HINT ${CGET_PREFIX}/lib)
find_library(BOOST_PYTHON_LIB boost_python27 HINT ${CGET_PREFIX}/lib)
find_library(PYTHON_LIB python2.7 HINT ${CGET_PREFIX}/lib)
find_library(ARMADILLO_LIB armadillo)
find_package(PythonLibs REQUIRED)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

set(SOURCE_FILES LDServer.cpp LDServer.h Raw.h Raw.cpp Segment.cpp Segment.h Cell.cpp Cell.h Types.h Morton.cpp Morton.h)
add_library(LDServer ${SOURCE_FILES})
target_link_libraries(LDServer ${SAVVY_LIB} ${HTS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB})

add_library(PyLDServer SHARED PyLDServer.cpp LDServer.cpp Raw.cpp Segment.cpp Segment.h Cell.cpp Cell.h Types.h Morton.cpp Morton.h)
set_target_properties(PyLDServer PROPERTIES SUFFIX ".so")
set_target_properties(PyLDServer PROPERTIES PREFIX "")
target_include_directories(PyLDServer PUBLIC ${PYTHON_INCLUDE_DIRS})
target_link_libraries(PyLDServer ${SAVVY_LIB} ${HTS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${BOOST_PYTHON_LIB} ${PYTHON_LIB})

install(TARGETS PyLDServer COMPONENT api EXPORT PyLDServer-config
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

add_subdirectory(test)