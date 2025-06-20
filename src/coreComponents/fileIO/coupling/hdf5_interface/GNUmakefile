HDF5_DIR      	 = /usr/tce/packages/hdf5/hdf5-parallel-1.8.18-gcc-4.9.3-mvapich2-2.2
HDF5INCFLAGS     = -I$(HDF5_DIR)/include
HDF5LIBFLAGS	 = -L$(HDF5_DIR)/lib

MPI_DIR 		 = /usr/tce/packages/mvapich2/mvapich2-2.2-gcc-4.9.3
HDFMPIINCFLAGS   = -I$(MPI_DIR)/include
HDFMPILIBFLAGS   = -L$(MPI_DIR)/lib -lhdf5 -lz


CXX = $(MPI_DIR)/bin/mpicxx -std=c++11
CXXFLAGS = -g
PPFLAGS = $(HDF5INCFLAGS) $(HDFMPIINCFLAGS)
LDFLAGS = 
LIBFLAGS = $(HDF5LIBFLAGS) $(HDFMPILIBFLAGS)

all: GEOSDriver.exe ChomboDriver.exe

%.exe: %.o GNUmakefile coupler.o
	$(CXX) $(CXXFLAGS) $(PPFLAGS)  $(LDFLAGS) $< coupler.o $(LIBFLAGS) -o $@

%.o: %.cpp GNUmakefile
	$(CXX) $(CXXFLAGS) -c $(PPFLAGS)  $< -o $@
	$(CXX) -MM $(PPFLAGS) $< > $*.d

vars:
	echo $(CXX)
	echo $(PPFLAGS)
	echo $(LIBFLAGS)
	echo $(LDFLAGS)

clean:
	rm -rf *.hdf5 *.exe *.o *.d *.dSYM

-include $(OBJS:.o=.d)
