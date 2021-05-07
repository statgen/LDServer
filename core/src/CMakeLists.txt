find_library(SAVVY_LIB savvy HINT ${CGET_PREFIX}/lib)
find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(TABIXPP_LIB tabixpp HINT ${CGET_PREFIX}/lib)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(ZSTD_LIB zstd HINT ${CGET_PREFIX}/lib)
find_library(HIREDIS_LIB hiredis HINT ${CGET_PREFIX}/lib)
find_library(ARMADILLO_LIB armadillo HINT ${CGET_PREFIX}/lib)
find_library(MKL_RT_LIB mkl_rt)
find_package(Python3 COMPONENTS Development REQUIRED)

set(BOOST_ROOT "" CACHE PATH ${CGET_PREFIX})
find_package(Boost REQUIRED COMPONENTS python36 iostreams)
find_package(Threads REQUIRED)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

if(NOT MKL_RT_LIB)
   set(MKL_RT_LIB "")
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(SOURCE_FILES LDServer.cpp LDServer.h ScoreServer.cpp ScoreServer.h Raw.h Raw.cpp Segment.cpp Segment.h ScoreSegment.cpp ScoreSegment.h Cell.cpp Cell.h Types.h Morton.cpp Morton.h Phenotypes.h Phenotypes.cpp Mask.cpp Mask.h SummaryStatisticsLoader.h SummaryStatisticsLoader.cpp ScoreCovarianceRunner.cpp ScoreCovarianceRunner.h)
add_library(LDServer ${SOURCE_FILES})
target_include_directories(LDServer PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(LDServer ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${TABIXPP_LIB} ${Boost_LIBRARIES} ${Python3_LIBRARIES})

add_library(pywrapper SHARED pywrapper.cpp LDServer.cpp ScoreServer.cpp Raw.cpp Segment.cpp Segment.h ScoreSegment.cpp ScoreSegment.h Cell.cpp Cell.h Types.h Morton.cpp Morton.h Phenotypes.h Phenotypes.cpp Mask.cpp Mask.h SummaryStatisticsLoader.h SummaryStatisticsLoader.cpp ScoreCovarianceRunner.cpp ScoreCovarianceRunner.h)
set_target_properties(pywrapper PROPERTIES SUFFIX ".so")
set_target_properties(pywrapper PROPERTIES PREFIX "")
target_include_directories(pywrapper PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(pywrapper ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${MKL_RT_LIB} ${TABIXPP_LIB} ${Boost_LIBRARIES} ${Python3_LIBRARIES})

add_executable(sandbox Sandbox.cpp ${SOURCE_FILES})
target_include_directories(sandbox PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(sandbox Threads::Threads ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${MKL_RT_LIB} ${TABIXPP_LIB} ${Boost_LIBRARIES} ${Python3_LIBRARIES})

install(TARGETS pywrapper DESTINATION ${PROJECT_SOURCE_DIR}/../rest/core)
