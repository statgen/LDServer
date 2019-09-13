cmake_minimum_required(VERSION 3.13)
project(libRmath VERSION 1.0.0)

add_custom_target(libRmath
  ALL COMMAND make CFLAGS="\"-I${CMAKE_PREFIX_PATH}/include ${CMAKE_CXX_FLAGS}\""
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  COMMENT "Building libRmath..."
)

install(DIRECTORY Rmath/ DESTINATION include)
install(FILES ${CMAKE_STATIC_LIBRARY_PREFIX}Rmath${CMAKE_STATIC_LIBRARY_SUFFIX} DESTINATION lib)
