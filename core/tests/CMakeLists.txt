find_library(GTEST_LIB gmock HINT ${CGET_PREFIX}/lib)
find_library(BOOST_SYSTEM_LIB boost_system HINT ${CGET_PREFIX}/lib)
find_library(BOOST_IOSTREAMS_LIB boost_iostreams HINT ${CGET_PREFIX}/lib)
find_library(MKL_RT_LIB mkl_rt)
find_package(Threads REQUIRED)

if(CGET_PREFIX)
    include_directories(${CGET_PREFIX}/include)
endif()

if(NOT MKL_RT_LIB)
   set(MKL_RT_LIB "")
endif()

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
add_executable(testAll LDServerTest.cpp Main_TestAll.cpp RareMetal.cpp)
target_link_libraries(testAll LDServer ${BOOST_SYSTEM_LIB} ${BOOST_IOSTREAMS_LIB} ${GTEST_LIB} ${MKL_RT_LIB} ${CMAKE_THREAD_LIBS_INIT})

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
        ../../data/chr21.test.tab
        ../../data/chr21.test.ped
        ../../data/chr21.test.dat
        ../../data/chr21.test.frq
        ../../data/chr22.test.vcf.gz
        ../../data/chr22.test.vcf.gz.tbi
        ../../data/chr22.test.bcf
        ../../data/chr22.test.bcf.csi
        ../../data/chr22.test.sav
        ../../data/chr22.test.sav.s1r
        ../../data/chr22.test.tab
        ../../data/chr22.test.bad_float.tab
        ../../data/chr22.test.missing_values.tab
        ../../data/chr22.test.ped
        ../../data/chr22.test.dat
        ../../data/chr22.test.frq
        ../../data/chr21.test.missing_genotypes_and_phenotypes.RAND_QT.singlevar.cov.txt
        ../../data/chr21.test.missing_genotypes_and_phenotypes.RAND_QT.singlevar.score.txt
        ../../data/chr21.test.missing_pheno.RAND_QT.singlevar.cov.txt
        ../../data/chr21.test.missing_pheno.RAND_QT.singlevar.score.txt
        ../../data/chr21.test.missing_values.ped
        ../../data/chr21.test.missing_values.tab
        ../../data/chr21.test.missing_values.vcf.gz
        ../../data/chr21.test.missing_values.vcf.gz.tbi
        ../../data/chr22.monomorphic_test.vcf.gz
        ../../data/chr22.monomorphic_test.vcf.gz.tbi
        ../../data/test_no_testable_variants.mask.tab.gz
        ../../data/test_no_testable_variants.tab
        ../../data/test_no_testable_variants.vcf.gz.tbi
        ../../data/test_no_testable_variants.mask.tab.gz.tbi
        ../../data/test_no_testable_variants.vcf.gz
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
        ../../data/mask.epacts.chr22.gencode-exons-AF01.tab.gz
        ../../data/mask.epacts.chr22.gencode-exons-AF01.tab.gz.tbi
        ../../data/mask.epacts.chr22.gencode-exons-AF05.tab.gz
        ../../data/mask.epacts.chr22.gencode-exons-AF05.tab.gz.tbi
        ../../data/chr21.test.RAND_QT.singlevar.cov.txt
        ../../data/chr21.test.RAND_QT.singlevar.score.txt
        DESTINATION test)
