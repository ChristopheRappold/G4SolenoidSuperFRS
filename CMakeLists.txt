#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(G4SolSimple)

include(CheckCXXCompilerFlag)
set(HAVE_ISYSTEM true)

set (CMAKE_CXX_STANDARD 17)

SET(VGM_INSTALL $ENV{VGM_INSTALL})
#SET(Geant4_DIR $ENV{GEANT4_DIR})
#SET(ROOT_DIR $ENV{ROOTSYS})
if(${VGM_INSTALL} MATCHES 4\\.2)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-4.2.0)
  message(STATUS "VGM dir 4.2:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 4\\.3)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-4.3.0)
  message(STATUS "VGM dir 4.3:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 4\\.4)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-4.4.0)
  message(STATUS "VGM dir 4.4:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 4\\.5)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-4.5.0)
  message(STATUS "VGM dir 4.5:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 4\\.8)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-4.8.0)
  message(STATUS "VGM dir 4.8:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 5\\.0)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-5.0.0)
  message(STATUS "VGM dir 5.0:" ${VGM_DIR})
elseif(${VGM_INSTALL} MATCHES 5\\.1)
  set(VGM_DIR ${VGM_INSTALL}/lib/VGM-5.1.0)
  message(STATUS "VGM dir 5.1:" ${VGM_DIR})
else()
  message(STATUS "no match for VGM version" ${VGM_INSTALL})
endif()

#set(CMAKE_MODULE_PATH 
#  ${CMAKE_CURRENT_SOURCE_DIR}/cmake
#  ${VGM_DIR}/Modules
#  ${CMAKE_MODULE_PATH}) 

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
message(STATUS "CMAKE_PREFIX_PATH " ${CMAKE_PREFIX_PATH})

find_package(ROOT REQUIRED)
if(ROOT_FOUND)
  add_definitions(-DCEXMC_USE_ROOT)
  include_directories(${ROOT_INCLUDE_DIRS})
  list(APPEND EXTRA_LIBRARIES ${ROOT_LIBRARIES})
  message(STATUS "ROOT Config :" ${ROOT_USE_FILE})
  message(STATUS "PATH: " $ENV{PATH})
endif()

find_package(VGM REQUIRED)
if (VGM_FOUND)
  message(STATUS "VGM :" ${VGM_INCLUDE_DIRS})
  add_definitions(-DG4SOLCONVERT)
  include_directories(${VGM_INCLUDE_DIRS})
endif()

find_package(Boost 1.54.0 REQUIRED)
if(Boost_FOUND)
  message(STATUS "BOOST :" ${Boost_INCLUDE_DIRS})
  include_directories(${Boost_INCLUDE_DIRS})
endif()
  
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
message(STATUS "Geant4 Setup:" ${CMAKE_CONFIG_PATH})

if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all qt)
  if(Geant4_FOUND)
    message(STATUS "Geant4 UI:" ${Geant4_INCLUDE_DIRS})
    message(STATUS "Geant4 def " ${Geant4_DEFINITIONS})
  endif()
else()
  find_package(Geant4 REQUIRED)
  if(Geant4_FOUND)
    message(STATUS "Geant4 :" ${Geant4_INCLUDE_DIRS})
  endif()
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
		    ${VGM_INCLUDE_DIRS}
		    )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(library_name vmc_${PROJECT_NAME})


add_executable(G4SolenoidSimple G4SolenoidSimple.cc ${sources} ${headers})
target_link_libraries(G4SolenoidSimple ${library_name} ${VGM_LIBRARIES} ${ROOT_LIBRARIES} ${Geant4_LIBRARIES})

ROOT_GENERATE_DICTIONARY(
  ${CMAKE_SHARED_LIBRARY_PREFIX}${library_name}_dict
#  no_rootmap
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_Hit.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_Event.hh
  LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/include/TG4Sol_LinkDef.h
  )
set(SourceLib ${PROJECT_SOURCE_DIR}/src/TG4Sol_Event.cc ${PROJECT_SOURCE_DIR}/src/TG4Sol_Hit.cc)
add_library(${library_name} SHARED ${SourceLib} ${CMAKE_SHARED_LIBRARY_PREFIX}${library_name}_dict.cxx ${headers})
target_link_libraries(${library_name} ${ROOT_LIBRARIES})

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
  runTemplate.mac
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

