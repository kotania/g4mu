#==============================================================================
#  CMakeLists.txt for building Geant4Py user modules
#==============================================================================
cmake_minimum_required(VERSION 3.16...3.21)

#==============================================================================
project(geometry)
#==============================================================================
# Geant4Py cmake path
set(geant4py_DIR /groups/icecube/kozynets/software/g4python/install/share/cmake)

# Geant4 installation path
set(GEANT4_INSTALL /groups/icecube/kozynets/software/geant4-v11.1.0-install
    CACHE STRING "Geant4 installation path")

# Pybind11 cmake files location
# Mac Homebrew
set(pybind11_DIR /cvmfs/icecube.opensciencegrid.org/py3-v4.2.1/RHEL_7_x86_64/share/cmake/pybind11
    CACHE STRING "Pybind11 search path")

# alternaive pybind11 path
#set(pybind11_DIR /opt/pybind11/share/cmake/pybind11
#    CACHE STRING "Pybind11 search path")

#------------------------------------------------------------------------------
# python module install path
set(MODULE ${CMAKE_PROJECT_NAME})
set(G4MODULES_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})

include(${geant4py_DIR}/geant4pyCommon.cmake)
#------------------------------------------------------------------------------
# add source
file(GLOB sources ${PROJECT_SOURCE_DIR}/${MODULE}/*.cc)

pybind11_add_module(${MODULE} ${sources})

target_include_directories(${MODULE} PRIVATE
  ${GEANT4_INCLUDE_DIR}
  ${Python3_INCLUDE_DIRS}
  ${pybind11_INCLUDE_DIR}
)

target_link_directories(${MODULE} PRIVATE ${GEANT4_LIBRARY_DIR})

target_link_libraries(${MODULE} PRIVATE
  ${GEANT4_LIBRARIES} ${Python3_LIBRARIES}
  global_cflags
)

set_target_properties(${MODULE} PROPERTIES PREFIX "")
set_target_properties(${MODULE} PROPERTIES OUTPUT_NAME ${MODULE})
set_target_properties(${MODULE} PROPERTIES SUFFIX ".so")

# install
install(TARGETS ${MODULE} LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
