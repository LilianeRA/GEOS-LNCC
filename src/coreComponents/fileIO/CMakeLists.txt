# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 TotalEnergies
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
#
#--------------------------------------------------------------------------------------------------

#[[
Package: fileIO

Contains:
 - I/O interfaces for the packing components and supported outputs and their python wrappers.
 - a coupler for data exchange with CHOMBO.
#]]

#
# Specify all headers
#
set( fileIO_headers
     LogLevelsInfo.hpp
     Outputs/BlueprintOutput.hpp
     Outputs/OutputBase.hpp
     Outputs/OutputManager.hpp
     Outputs/OutputUtilities.hpp
     Outputs/PythonOutput.hpp
     Outputs/RestartOutput.hpp
     Outputs/TimeHistoryOutput.hpp
     timeHistory/HDFFile.hpp
     timeHistory/HistoryCollectionBase.hpp
     timeHistory/BufferedHistoryIO.hpp
     timeHistory/PackCollection.hpp
     timeHistory/HDFHistoryIO.hpp
     timeHistory/HistoryCollection.hpp )

#
# Specify all sources
#
set( fileIO_sources
     Outputs/BlueprintOutput.cpp
     Outputs/OutputBase.cpp
     Outputs/OutputManager.cpp
     Outputs/OutputUtilities.cpp
     Outputs/PythonOutput.cpp
     Outputs/RestartOutput.cpp
     Outputs/TimeHistoryOutput.cpp
     timeHistory/HDFFile.cpp
     timeHistory/HistoryCollectionBase.cpp
     timeHistory/PackCollection.cpp
     timeHistory/HDFHistoryIO.cpp )

set( dependencyList ${parallelDeps} events mesh HDF5::HDF5 )
if( ENABLE_PYGEOSX )
  list( APPEND fileIO_headers
        python/PyHistoryCollectionType.hpp
        python/PyHistoryOutputType.hpp
        python/PyVTKOutputType.hpp )
  list( APPEND fileIO_sources
        python/PyHistoryCollection.cpp
        python/PyHistoryOutput.cpp
        python/PyVTKOutput.cpp )
  list( APPEND dependencyList Python3::Python pylvarray )
endif()


if( ENABLE_MPI )
  add_subdirectory( coupling/hdf5_interface )
  list( APPEND dependencyList mpi hdf5_interface HDF5::HDF5 )
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp Outputs/ChomboIO.hpp )
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp Outputs/ChomboIO.cpp )
endif()

if( ENABLE_SILO )
  list( APPEND dependencyList silo HDF5::HDF5)
  list( APPEND fileIO_headers
               silo/SiloFile.hpp
               Outputs/SiloOutput.hpp )
  list( APPEND fileIO_sources
               silo/SiloFile.cpp
               Outputs/SiloOutput.cpp )
endif( )

if( ENABLE_VTK )
  list( APPEND fileIO_headers
        vtk/VTKPVDWriter.hpp
        vtk/VTKVTMWriter.hpp
        vtk/VTKPolyDataWriterInterface.hpp
        Outputs/VTKOutput.hpp
        )
  list( APPEND fileIO_sources
        vtk/VTKPVDWriter.cpp
        vtk/VTKVTMWriter.cpp
        vtk/VTKPolyDataWriterInterface.cpp
        Outputs/VTKOutput.cpp
        )
    list( APPEND dependencyList VTK::IOLegacy VTK::IOXML )
endif()

if( ENABLE_CUDA AND ENABLE_CUDA_NVTOOLSEXT )
  list( APPEND dependencyList CUDA::nvToolsExt )
endif()

geos_decorate_link_dependencies( LIST decoratedDependencies
                                 DEPENDENCIES ${dependencyList} )

blt_add_library( NAME       fileIO
                 SOURCES    ${fileIO_sources}
                 HEADERS    ${fileIO_headers}
                 DEPENDS_ON ${decoratedDependencies}
                 OBJECT     ${GEOS_BUILD_OBJ_LIBS}
                 SHARED     ${GEOS_BUILD_SHARED_LIBS}
               )


target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

install( TARGETS fileIO LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )



