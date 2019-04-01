cmake_minimum_required(VERSION 3.10)
project(FWH CXX)

include_directories("${PROJECT_BINARY_DIR}/include")

if(NOT CMAKE_BUILD_TYPE) # force set build type in cache to debug
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# force set gui setting
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -D_DEBUG" CACHE STRING "Flags used by the linker during debug builds" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall" CACHE STRING "Flags used by the linker during release builds" FORCE)

#HDF5 library
set(HDF5_INCLUDE "/opt" CACHE PATH "path to hdf5 include")
set(HDF5_LD "/opt" CACHE PATH "path to hdf5 library")
include_directories(${HDF5_INCLUDE})
set(CXX_LD ${CXX_LD} ${HDF5_LD})
set(CXX_LIB ${CXX_LIB} hdf5)


#source file list
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRCLIST)

# use C++14
add_definitions(-std=c++14)

#build
LINK_DIRECTORIES(${CXX_LD})
add_executable(FWH ${SRCLIST})
target_link_libraries(FWH ${CXX_LIB})
