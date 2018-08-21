cmake_minimum_required(VERSION 3.8)
project(LDServer VERSION 0.0.1)
include(GNUInstallDirs)
add_compile_options(-Wall -Wpedantic)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)

add_subdirectory(src)
add_subdirectory(tests)