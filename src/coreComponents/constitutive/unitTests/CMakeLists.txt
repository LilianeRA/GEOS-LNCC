# Specify list of tests
set( gtest_geosx_tests
     testBrineSalinity.cpp
     testCompositionalDensity.cpp
     testCompositionalPhillipsBrineDensity.cpp
     testCompositionalPhillipsBrineViscosity.cpp
     testCompositionalProperties.cpp
     testCompositionalFluidComponentType.cpp
     testCubicEOS.cpp
     testDamageUtilities.cpp
     testDruckerPrager.cpp
     testElasticIsotropic.cpp
     testKValueInitialization.cpp
     testImmiscibleWaterFlashModel.cpp
     testImmiscibleWaterProperties.cpp
     testLohrenzBrayClarkViscosity.cpp
     testModifiedCamClay.cpp
     testMultiFluidBlackOil.cpp
     testMultiFluidCO2Brine.cpp
     testMultiFluidTwoPhaseCompositionalMultiphase.cpp
     testMultiFluidDeadOil.cpp
     testMultiFluidSelector.cpp
     testNegativeTwoPhaseFlash.cpp
     testNegativeTwoPhaseFlash9Comp.cpp
     testParticleFluidEnums.cpp
     testPressureTemperatureCoordinates.cpp
     testPropertyConversions.cpp
     testSoreideWhitsonEOSPhaseModel.cpp
     testSoreideWhitsonFlash.cpp
     testStabilityTest2Comp.cpp
     testStabilityTest9Comp.cpp
     testRachfordRice.cpp
     testTwoPhaseImmiscibleFluid.cpp )

set( dependencyList gtest blas lapack constitutive ${parallelDeps} )

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

# Add gtest C++ based tests
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${dependencyList} )

    geos_add_test( NAME ${test_name}
                   COMMAND ${test_name} )

endforeach()
