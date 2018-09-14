find_library(GTEST_LIB gmock HINT ${CGET_PREFIX}/lib)

find_package(Threads REQUIRED)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

add_executable(testAll LDServerTest.cpp Main_TestAll.cpp)
target_link_libraries(testAll LDServer ${GTEST_LIB} ${CMAKE_THREAD_LIBS_INIT})

install(TARGETS testAll COMPONENT cli RUNTIME DESTINATION test OPTIONAL)
install(FILES
        ../../data/EUR.samples.txt
        ../../data/AFR.samples.txt
        ../../data/chr21.test.vcf.gz
        ../../data/chr21.test.vcf.gz.tbi
        ../../data/chr21.test.bcf
        ../../data/chr21.test.bcf.csi
        ../../data/chr21.test.sav
        ../../data/chr21.test.sav.s1r
        ../../data/chr22.test.vcf.gz
        ../../data/chr22.test.vcf.gz.tbi
        ../../data/chr22.test.bcf
        ../../data/chr22.test.bcf.csi
        ../../data/chr22.test.sav
        ../../data/chr22.test.sav.s1r
        ../../data/chrX.test.sav
        ../../data/chrX.test.sav.s1r
        ../../data/region_ld_22_51241101_51241385.hap.ld
        ../../data/region_ld_22_50544251_50549251.hap.ld
        ../../data/region_ld_22_51241101_51241385.AFR.hap.ld
        ../../data/variant_ld_22_51241101_vs_51241101_51241385.hap.ld
        ../../data/variant_ld_22_51241386_vs_51241101_51241385.hap.ld
        ../../data/variant_ld_22_51241309_vs_51241101_51244237.hap.ld
        ../../data/variant_ld_22_50546666_vs_50544251_50549251.hap.ld
        ../../data/region_ld_X_60100_60150.hap.ld
        DESTINATION test)