set( subdirs
     math
     common
     codingUtilities
     dataRepository
     denseLinearAlgebra
     functions
     constitutive
     schema
     finiteElement
     mesh
     fieldSpecification
     finiteVolume
     linearAlgebra
     discretizationMethods
     events
     constitutiveDrivers
     fileIO
     physicsSolvers
     mainInterface )

unset( parallelDeps )

if ( ENABLE_OPENMP )
  list( APPEND parallelDeps openmp )
endif()

if ( ENABLE_CUDA )
  list( APPEND parallelDeps cuda )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND parallelDeps CUDA::nvToolsExt )
endif()

if( ENABLE_HIP )
  list( APPEND parallelDeps blt::hip )
endif()

if( ENABLE_MPI )
  list( APPEND parallelDeps mpi )
endif()

if ( GEOS_BUILD_OBJ_LIBS )
  set( LVARRAY_BUILD_OBJ_LIBS TRUE CACHE BOOL "" FORCE )
endif()

add_subdirectory( LvArray )

if( ENABLE_PYGEOSX )
  list( APPEND coreLibs pylvarray )
endif()

foreach( lib ${subdirs} )
  if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
    list( APPEND coreLibs ${lib} )
  endif()
endforeach()

foreach( lib ${subdirs} )
  add_subdirectory( ${lib} )
endforeach()

add_library(geosx_core INTERFACE)

if( GEOS_BUILD_SHARED_LIBS)
  target_link_libraries(geosx_core INTERFACE mainInterface )
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_link_options(geosx_core INTERFACE -Wl,--no-as-needed)
  endif()
else()
  target_link_libraries(geosx_core INTERFACE $<LINK_LIBRARY:WHOLE_ARCHIVE,mainInterface> )
endif()

geosx_add_code_checks( PREFIX coreComponents
                       EXCLUDES cmake constitutive/PVTPackage )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check 2>&1 >/dev/null"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

if( GEOS_ENABLE_TESTS )
  add_subdirectory( unitTests )
endif( )

