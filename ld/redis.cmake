cmake_minimum_required(VERSION 3.8)
project(redis VERSION 5.0.4)

execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
install(PROGRAMS src/redis-server src/redis-cli DESTINATION bin)