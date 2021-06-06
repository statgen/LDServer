find_library(SAVVY_LIB savvy HINT ${CGET_PREFIX}/lib)
find_library(HTS_LIB hts HINT ${CGET_PREFIX}/lib)
find_library(TABIXPP_LIB tabixpp HINT ${CGET_PREFIX}/lib)
find_library(Z_LIB NAMES libz.a z HINT ${CGET_PREFIX}/lib)
find_library(ZSTD_LIB zstd HINT ${CGET_PREFIX}/lib)
find_library(HIREDIS_LIB hiredis HINT ${CGET_PREFIX}/lib)
find_library(ARMADILLO_LIB armadillo HINT ${CGET_PREFIX}/lib)
find_library(MKL_RT_LIB mkl_rt)
find_library(ARROW_LIB arrow HINT ${CGET_PREFIX/lib})
find_library(PARQUET_LIB parquet HINT ${CGET_PREFIX/lib})
find_package(Python3 3.8 COMPONENTS Development REQUIRED)

set(BOOST_ROOT "" CACHE PATH ${CGET_PREFIX})
find_package(Boost REQUIRED COMPONENTS python38 iostreams)
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
target_link_libraries(LDServer ${PARQUET_LIB} ${ARROW_LIB} ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${Boost_LIBRARIES} ${TABIXPP_LIB})

add_library(pywrapper SHARED pywrapper.cpp LDServer.cpp ScoreServer.cpp Raw.cpp Segment.cpp Segment.h ScoreSegment.cpp ScoreSegment.h Cell.cpp Cell.h Types.h Morton.cpp Morton.h Phenotypes.h Phenotypes.cpp Mask.cpp Mask.h SummaryStatisticsLoader.h SummaryStatisticsLoader.cpp ScoreCovarianceRunner.cpp ScoreCovarianceRunner.h)
set_target_properties(pywrapper PROPERTIES SUFFIX ".so")
set_target_properties(pywrapper PROPERTIES PREFIX "")
target_include_directories(pywrapper PUBLIC ${Python3_INCLUDE_DIRS})
target_link_libraries(pywrapper ${PARQUET_LIB} ${ARROW_LIB} ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${MKL_RT_LIB} ${TABIXPP_LIB} ${Boost_LIBRARIES} ${Python3_LIBRARIES})

add_executable(sandbox Sandbox.cpp ${SOURCE_FILES})
target_link_libraries(sandbox Threads::Threads ${PARQUET_LIB} ${ARROW_LIB} ${SAVVY_LIB} ${HTS_LIB} ${HIREDIS_LIB} ${Z_LIB} ${ZSTD_LIB} ${ARMADILLO_LIB} ${MKL_RT_LIB} ${TABIXPP_LIB} ${Boost_LIBRARIES} ${Python3_LIBRARIES})

install(TARGETS pywrapper DESTINATION ${PROJECT_SOURCE_DIR}/../rest/core)
