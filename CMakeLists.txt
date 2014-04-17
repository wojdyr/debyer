cmake_minimum_required(VERSION 2.8.5 FATAL_ERROR)
project("Debyer")

add_definitions(-DVERSION="0.3")

find_package(ZLIB)
find_package(BZip2)
find_package(MPI)

option(ENABLE_SINGLE "Use single precision" off)
option(ENABLE_ZLIB   "Use zlib" on)
option(ENABLE_BZLIB  "Use bzip2" on)
option(ENABLE_MPI    "Enable MPI" off)

# single precision?
if(ENABLE_SINGLE)
  add_definitions(-DUSE_SINGLE)
endif()

# zlib
if(ENABLE_ZLIB AND ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  set(LINKLIBS ${LINKLIBS} ${ZLIB_LIBRARIES})
  add_definitions(-DHAVE_ZLIB)
endif()

# bzip2
if(ENABLE_BZLIB AND BZIP2_FOUND)
  include_directories(${BZIP2_INCLUDE_DIR})
  set(LINKLIBS ${LINKLIBS} ${BZIP2_LIBRARIES})
  add_definitions(-DHAVE_BZLIB)
endif()

# MPI
if(ENABLE_MPI AND MPI_CXX_FOUND)
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
  if (APPLE)
    set(CMAKE_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS})
  else()
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  endif()
  set(LINKLIBS ${LINKLIBS} ${MPI_CXX_LIBRARIES})
  include_directories(${MPI_CXX_INCLUDE_PATH})
  add_definitions(-DENABLE_MPI)
endif()

add_subdirectory(debyer)

