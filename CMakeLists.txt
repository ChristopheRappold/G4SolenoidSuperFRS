#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(G4SolSimple)

include(CheckCXXCompilerFlag)
set(HAVE_ISYSTEM true)

set (CMAKE_CXX_STANDARD 14)

SET(Geant4_DIR $ENV{GEANT4_DIR})
SET(ROOT_DIR $ENV{ROOTSYS})

set(CMAKE_MODULE_PATH 
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
  ${CMAKE_MODULE_PATH}) 


find_package(ROOT REQUIRED)
message(STATUS "ROOT Config :" ${ROOT_USE_FILE})
message(STATUS "PATH: " $ENV{PATH})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
		    )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(library_name vmc_${PROJECT_NAME})


add_executable(G4SolenoidSimple G4SolenoidSimple.cc ${sources} ${headers})
target_link_libraries(G4SolenoidSimple ${library_name} ${ROOT_LIBRARIES} ${Geant4_LIBRARIES})

ROOT_GENERATE_DICTIONARY(
  ${CMAKE_SHARED_LIBRARY_PREFIX}${library_name}
  no_rootmap
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_Hit.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_Event.hh
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_LinkDef.h
  )
set(SourceLib ${PROJECT_SOURCE_DIR}/src/TG4Sol_Event.cc ${PROJECT_SOURCE_DIR}/src/TG4Sol_Hit.cc)
add_library(${library_name} ${SourceLib} ${CMAKE_SHARED_LIBRARY_PREFIX}${library_name}_dict.cxx ${headers})
target_link_libraries(${library_name})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build G4SolenoidSimple. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(G4SOL_SCRIPTS
  G4SolenoidSimple.in 
#  G4SolenoidSimple.out 
  testconfig.par
  icons.mac
  gui.mac
  run.png
  init.mac 
  init_vis.mac 
  run1.mac 
  run2.mac 
  vis.mac
  )

foreach(_script ${G4SOL_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(G4SolSimple DEPENDS G4SolenoidSimple)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS G4SolenoidSimple DESTINATION bin)

