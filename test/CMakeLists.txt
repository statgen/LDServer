find_library(GTEST_LIB gtest HINT ${CGET_PREFIX}/lib)
find_package(Threads REQUIRED)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

add_executable(testAll LDServerTest.cpp Main_TestAll.cpp)
target_link_libraries(testAll LDServer ${GTEST_LIB} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS testAll COMPONENT cli RUNTIME DESTINATION test OPTIONAL)
install(FILES
        chr21.test.vcf.gz
        chr21.test.vcf.gz.tbi
        chr21.test.bcf
        chr21.test.bcf.csi
        chr21.test.sav
        chr21.test.sav.s1r
        chr22.test.vcf.gz
        chr22.test.vcf.gz.tbi
        chr22.test.bcf
        chr22.test.bcf.csi
        chr22.test.sav
        chr22.test.sav.s1r
        chrX.test.sav
        chrX.test.sav.s1r
        region_ld_22_51241101_51241385.hap.ld
        variant_ld_22_51241101_vs_51241101_51241385.hap.ld
        region_ld_X_60100_60150.hap.ld
        DESTINATION test)