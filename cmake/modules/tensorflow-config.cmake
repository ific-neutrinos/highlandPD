# Locates the tensorFlow library and include directories.

IF(EXISTS /cvmfs/larsoft.opensciencegrid.org/)
          set(TENSORFLOW_LIBRARIES /cvmfs/larsoft.opensciencegrid.org/products/tensorflow/v1_12_0c/Linux64bit+2.6-2.12-e19-py2-prof/lib/libtensorflow-core.so)
          set(TENSORFLOW_INCLUDE_DIRS /cvmfs/larsoft.opensciencegrid.org/products/tensorflow/v1_12_0c/Linux64bit+2.6-2.12-e19-py2-prof/include/ /cvmfs/larsoft.opensciencegrid.org/products/tensorflow/v1_12_0c/Linux64bit+2.6-2.12-e19-py2-prof/tensorflow-1.12.0/tensorflow/contrib/makefile/downloads/absl/)
	  set(COMPILETF 1)
else()
	  set(COMPILETF 0)
endif()
