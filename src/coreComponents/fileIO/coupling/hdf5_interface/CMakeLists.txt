#
# Specify all headers
#
set(hdf5_interface_headers
  coupler.hpp )

#
# Specify all sources
#
set(hdf5_interface_sources
  coupler.cpp )

set( dependencyList HDF5::HDF5 mpi )

blt_add_library( NAME                  hdf5_interface
                 SOURCES               ${hdf5_interface_sources}
                 HEADERS               ${hdf5_interface_headers}
                 DEPENDS_ON            ${dependencyList}
               )
               
target_include_directories( hdf5_interface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

install( TARGETS hdf5_interface LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib )
